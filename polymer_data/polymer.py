#!/usr/bin/env python3
# Linear Polymer Generator for LAMMPS
# Creates a linear polymer chain that folds into a rectangular pattern when reaching box boundaries
# Modified to generate multiple polymer systems with different numbers of atoms

import numpy as np
import argparse
import os

def generate_polymer(num_atoms, box_size, bond_length=1.0):
    """
    Generate a linear polymer chain that folds into a rectangular pattern
    
    Parameters:
    -----------
    num_atoms : int
        Number of atoms in the polymer
    box_size : list
        Box dimensions [xlo, xhi, ylo, yhi, zlo, zhi]
    bond_length : float
        Bond length between consecutive atoms
    
    Returns:
    --------
    positions : ndarray
        Array of atom positions
    bonds : ndarray
        Array of bonds connecting atoms
    """
    # Extract box dimensions
    xlo, xhi, ylo, yhi, zlo, zhi = box_size
    
    # Width and height of the box
    width = xhi - xlo - 2*bond_length  # Leave some margin
    height = yhi - ylo - 2*bond_length  # Leave some margin
    
    # Initialize positions array
    positions = np.zeros((num_atoms, 3))
    
    # Initial position (with some margin from the boundary)
    x, y, z = xlo + bond_length, ylo + bond_length, (zlo + zhi) / 2
    
    # Direction flags (1: positive, -1: negative)
    x_direction = 1
    
    # Generate positions for each atom
    for i in range(num_atoms):
        positions[i] = [x, y, z]
        
        # Advance position along x-axis
        x += x_direction * bond_length
        
        # Check if we hit the box boundary
        if x > xhi - bond_length or x < xlo + bond_length:
            # Change direction and move up
            x_direction *= -1
            x += x_direction * bond_length
            y += bond_length
            
            # Check if we hit the top boundary
            if y > yhi - bond_length:
                # Move down to the bottom
                y = ylo + bond_length
    
    # Generate bonds connecting consecutive atoms
    bonds = np.zeros((num_atoms - 1, 2), dtype=int)
    for i in range(num_atoms - 1):
        bonds[i] = [i + 1, i + 2]  # LAMMPS uses 1-based indexing
    
    return positions, bonds

def write_data_file(filename, num_atoms, box_size, positions, bonds):
    """
    Write LAMMPS data file
    
    Parameters:
    -----------
    filename : str
        Output filename
    num_atoms : int
        Number of atoms
    box_size : list
        Box dimensions [xlo, xhi, ylo, yhi, zlo, zhi]
    positions : ndarray
        Array of atom positions
    bonds : ndarray
        Array of bonds
    """
    with open(filename, 'w') as f:
        # Header
        f.write("LAMMPS Data File for Linear Polymer Chain\n\n")
        
        # System dimensions
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_atoms - 1} bonds\n")
        f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        
        # Atom and bond types
        f.write("1 atom types\n")
        f.write("1 bond types\n\n")
        
        # Box dimensions
        xlo, xhi, ylo, yhi, zlo, zhi = box_size
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")
        
        # Masses
        f.write("Masses\n\n")
        f.write("1 1.0  # Type A\n\n")
        
        # Atoms section - molecular format
        f.write("Atoms # molecular\n\n")
        for i in range(num_atoms):
            # id mol-id type x y z
            mol_id = 1  # Single molecule
            atom_type = 1  # Type A
            x, y, z = positions[i]
            f.write(f"{i+1} {mol_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\n")
        
        # Bonds section
        f.write("Bonds\n\n")
        for i in range(len(bonds)):
            # id type atom1 atom2
            bond_type = 1
            atom1, atom2 = bonds[i]
            f.write(f"{i+1} {bond_type} {atom1} {atom2}\n")

def main():
    parser = argparse.ArgumentParser(description='Generate a LAMMPS data file for a linear polymer')
    parser.add_argument('-n', '--num_atoms', type=int, nargs='+', default=[100], 
                        help='Number of atoms in the polymer chain(s). Multiple values can be specified.')
    parser.add_argument('-o', '--output', type=str, default='polymer.data', 
                        help='Output data file name (ignored if multiple -n values are provided)')
    parser.add_argument('-b', '--box_size', type=float, nargs=6, default=[0, 32.0, 0, 32.0, 0, 32.0],
                       help='Box dimensions: xlo xhi ylo yhi zlo zhi')
    parser.add_argument('-l', '--bond_length', type=float, default=1.0, 
                        help='Bond length between consecutive atoms')
    
    args = parser.parse_args()
    
    # Process each requested system size
    for n in args.num_atoms:
        # Generate the output filename
        if len(args.num_atoms) > 1 or args.output == 'polymer.data':
            output_file = f"polymer_system_{n}.data"
        else:
            output_file = args.output
            
        # Generate polymer positions and bonds
        positions, bonds = generate_polymer(n, args.box_size, args.bond_length)
        
        # Write LAMMPS data file
        write_data_file(output_file, n, args.box_size, positions, bonds)
        
        print(f"Generated LAMMPS data file '{output_file}' with {n} atoms and {n - 1} bonds")

if __name__ == "__main__":
    main()

# python3 polymer.py -n 50 100 150 200 250 300 