// Computes MSD and RG for a LAMMPS trajectory - through COM

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

struct Atom {
    int id, type;
    double x, y, z;
};

struct Frame {
    int timestep;
    vector<Atom> atoms;
};

// Calculate center of mass for a frame
vector<double> calculateCOM(const Frame& frame) {
    vector<double> com = {0.0, 0.0, 0.0};
    int num_atoms = frame.atoms.size();
    
    for (const auto& atom : frame.atoms) {
        com[0] += atom.x;
        com[1] += atom.y;
        com[2] += atom.z;
    }
    
    com[0] /= num_atoms;
    com[1] /= num_atoms;
    com[2] /= num_atoms;
    
    return com;
}

// Calculate radius of gyration for a frame
double calculateRg(const Frame& frame, const vector<double>& com) {
    double rg_sq = 0.0;
    int num_atoms = frame.atoms.size();
    
    for (const auto& atom : frame.atoms) {
        double dx = atom.x - com[0];
        double dy = atom.y - com[1];
        double dz = atom.z - com[2];
        rg_sq += (dx*dx + dy*dy + dz*dz);
    }
    
    rg_sq /= num_atoms;
    return sqrt(rg_sq);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <trajectory_file> <output_prefix>" << endl;
        return 1;
    }
    
    string trajectory_file = argv[1];
    string output_prefix = argv[2];
    
    ifstream infile(trajectory_file);
    if (!infile) {
        cerr << "Error: Unable to open trajectory file " << trajectory_file << endl;
        return 1;
    }
    
    ofstream msd_outfile(output_prefix + "_msd.dat");
    ofstream rg_outfile(output_prefix + "_rg.dat");
    
    if (!msd_outfile || !rg_outfile) {
        cerr << "Error: Unable to create output files" << endl;
        return 1;
    }
    
    msd_outfile << "# Timestep MSD" << endl;
    rg_outfile << "# Timestep Rg" << endl;
    
    string line;
    vector<Frame> frames;
    Frame current_frame;
    int num_atoms = 0;
    int line_count = 0;
    bool reading_atoms = false;
    int atom_count = 0;
    
    while (getline(infile, line)) {
        line_count++;
        
        if (line.find("ITEM: TIMESTEP") != string::npos) {
            if (!current_frame.atoms.empty()) {
                frames.push_back(current_frame);
            }
            
            current_frame = Frame();
            reading_atoms = false;
            atom_count = 0;
            getline(infile, line);
            current_frame.timestep = stoi(line);
        } 
        else if (line.find("ITEM: NUMBER OF ATOMS") != string::npos) {
            getline(infile, line);
            num_atoms = stoi(line);
            current_frame.atoms.reserve(num_atoms);
        }
        else if (line.find("ITEM: ATOMS") != string::npos) {
            reading_atoms = true;
            continue;
        }
        
        if (reading_atoms) {
            istringstream iss(line);
            Atom atom;
            if (!(iss >> atom.id >> atom.type >> atom.x >> atom.y >> atom.z)) {
                cerr << "Error parsing atom line: " << line_count << endl;
                continue;
            }
            current_frame.atoms.push_back(atom);
            atom_count++;
            
            if (atom_count == num_atoms) {
                reading_atoms = false;
            }
        }
    }
    
    // Add the last frame if it exists
    if (!current_frame.atoms.empty()) {
        frames.push_back(current_frame);
    }
    
    if (frames.empty()) {
        cerr << "Error: No frames were read from the trajectory file" << endl;
        return 1;
    }
    
    cout << "Read " << frames.size() << " frames from trajectory file" << endl;
    
    // Calculate initial center of mass
    vector<double> initial_com = calculateCOM(frames.front());
    
    // Calculate MSD and Rg for each frame
    for (const auto& frame : frames) {
        vector<double> current_com = calculateCOM(frame);
        
        double dx = current_com[0] - initial_com[0];
        double dy = current_com[1] - initial_com[1];
        double dz = current_com[2] - initial_com[2];
        
        double msd = dx*dx + dy*dy + dz*dz;
        double rg = calculateRg(frame, current_com);
        
        msd_outfile << frame.timestep << " " << msd << endl;
        rg_outfile << frame.timestep << " " << rg << endl;
    }
    
    msd_outfile.close();
    rg_outfile.close();
    
    cout << "Analysis complete. Results written to " << output_prefix << "_msd.dat and " 
              << output_prefix << "_rg.dat" << endl;
    
    return 0;
}