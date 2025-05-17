import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(x, a, b):
    return a * (x ** b)

def linear_fit(x, a, b):
    return a * x + b

def read_data_file(filename):
    data = np.loadtxt(filename, skiprows=1)
    return data[:, 0], data[:, 1]  # time, value

def process_msd_data(base_dir, n_values, seeds):
    """Process MSD data for all N values and seeds"""
    results = {}
    
    for n in n_values:
        n_results = []
        for seed in seeds:
            filepath = os.path.join(base_dir, f"N{n}_seed{seed}", f"N{n}_seed{seed}_msd.dat")
            print(f"Looking for MSD file: {filepath}")
            if os.path.exists(filepath):
                print(f"Found MSD file: {filepath}")
                time, msd = read_data_file(filepath)
                n_results.append((time, msd))
        
        if n_results:
            # Use the shortest time series length
            min_length = min(len(t) for t, _ in n_results)
            time = n_results[0][0][:min_length]
            msd_values = np.array([m[:min_length] for _, m in n_results])
            avg_msd = np.mean(msd_values, axis=0)
            std_msd = np.std(msd_values, axis=0)
            results[n] = (time, avg_msd, std_msd)
    
    return results

def process_rg_data(base_dir, n_values, seeds):
    """Process radius of gyration data for all N values and seeds"""
    results = {}
    
    for n in n_values:
        rg_values = []
        for seed in seeds:
            filepath = os.path.join(base_dir, f"N{n}_seed{seed}", f"N{n}_seed{seed}_rg.dat")
            print(f"Looking for Rg file: {filepath}")
            if os.path.exists(filepath):
                print(f"Found Rg file: {filepath}")
                time, rg = read_data_file(filepath)
                # Use the last 25% of data points for steady-state Rg
                steady_state_idx = int(len(rg))
                steady_rg = np.mean(rg[steady_state_idx])
                rg_values.append(steady_rg)
        
        if rg_values:
            avg_rg = np.mean(rg_values)
            std_rg = np.std(rg_values)
            results[n] = (avg_rg, std_rg)
    
    return results

def calculate_diffusion_coefficients(msd_results):
    """Calculate diffusion coefficients from MSD data"""
    diff_coeff = {}
    
    for n, (time, msd, std_msd) in msd_results.items():
        # Use linear region (skip initial transient)
        start_idx = 0
        end_idx = int(len(time))
        
        # MSD = 6*D*t for 3D diffusion
        slope, _ = np.polyfit(time[start_idx:end_idx], msd[start_idx:end_idx], 1)
        D = slope / 6.0
        
        # Estimate error from standard deviation
        slopes = []
        for offset in [-0.1, 0, 0.1]:  # Use variations around the mean
            if offset == 0:
                continue
            slope, _ = np.polyfit(time[start_idx:end_idx], 
                                 msd[start_idx:end_idx] + offset * std_msd[start_idx:end_idx], 1)
            slopes.append(slope / 6.0)
        
        D_err = np.std(slopes + [D])
        diff_coeff[n] = (D, D_err)
    
    return diff_coeff

def plot_msd_vs_time(msd_results, output_dir):
    """Plot MSD vs time for all N values"""
    plt.figure(figsize=(10, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(msd_results)))
    
    for i, (n, (time, msd, std_msd)) in enumerate(sorted(msd_results.items())):
        plt.plot(time, msd, label=f'N = {n}', color=colors[i])
        plt.fill_between(time, msd - std_msd, msd + std_msd, alpha=0.3, color=colors[i])
    
    plt.xlabel('Time')
    plt.ylabel('MSD')
    plt.title('Mean Square Displacement vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'msd_vs_time.png'), dpi=300)
    
    plt.figure(figsize=(10, 6))
    for i, (n, (time, msd, std_msd)) in enumerate(sorted(msd_results.items())):
        plt.loglog(time, msd, label=f'N = {n}', color=colors[i])
    
    plt.xlabel('Time (log scale)')
    plt.ylabel('MSD (log scale)')
    plt.title('Log-Log Plot of MSD vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    plt.savefig(os.path.join(output_dir, 'log_msd_vs_time.png'), dpi=300)

def plot_diffusion_vs_n(diff_coeff, output_dir):
    """Plot diffusion coefficient vs N and fit power law"""
    n_values = np.array(list(diff_coeff.keys()))
    d_values = np.array([d for d, _ in diff_coeff.values()])
    d_errs = np.array([err for _, err in diff_coeff.values()])
    
    plt.figure(figsize=(8, 6))
    plt.errorbar(n_values, d_values, yerr=d_errs, fmt='o', capsize=5)
    
    # Fit power law: D ~ N^alpha
    try:
        popt, pcov = curve_fit(power_law, n_values, d_values, sigma=d_errs, absolute_sigma=True)
        a, alpha = popt
        alpha_err = np.sqrt(np.diag(pcov))[1]
        
        x_fit = np.linspace(min(n_values) * 0.9, max(n_values) * 1.1, 100)
        y_fit = power_law(x_fit, a, alpha)
        
        plt.plot(x_fit, y_fit, 'r-', label=f'Fit: D ~ N^{alpha:.3f}±{alpha_err:.3f}')
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('N (chain length)')
        plt.ylabel('Diffusion Coefficient D')
        plt.title('Diffusion Coefficient vs Chain Length')
        plt.legend()
        plt.grid(True, alpha=0.3, which='both')
        plt.savefig(os.path.join(output_dir, 'diffusion_vs_n.png'), dpi=300)
        
        print(f"Diffusion coefficient scaling: D ~ N^{alpha:.3f} ± {alpha_err:.3f}")
    except Exception as e:
        print(f"Error fitting diffusion vs N: {e}")

def plot_rg_vs_n(rg_results, output_dir):
    """Plot radius of gyration vs N and fit power law"""
    n_values = np.array(list(rg_results.keys()))
    rg_values = np.array([rg for rg, _ in rg_results.values()])
    rg_errs = np.array([err for _, err in rg_results.values()])
    
    plt.figure(figsize=(8, 6))
    plt.errorbar(n_values, rg_values, yerr=rg_errs, fmt='o', capsize=5)
    
    # Fit power law: Rg ~ N^alpha
    try:
        popt, pcov = curve_fit(power_law, n_values, rg_values, sigma=rg_errs, absolute_sigma=True)
        a, alpha = popt
        alpha_err = np.sqrt(np.diag(pcov))[1]
        
        x_fit = np.linspace(min(n_values) * 0.9, max(n_values) * 1.1, 100)
        y_fit = power_law(x_fit, a, alpha)
        
        plt.plot(x_fit, y_fit, 'r-', label=f'Fit: Rg ~ N^{alpha:.3f}±{alpha_err:.3f}')
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('N (chain length)')
        plt.ylabel('Radius of Gyration Rg')
        plt.title('Radius of Gyration vs Chain Length')
        plt.legend()
        plt.grid(True, alpha=0.3, which='both')
        plt.savefig(os.path.join(output_dir, 'rg_vs_n.png'), dpi=300)
        
        print(f"Radius of gyration scaling: Rg ~ N^{alpha:.3f} ± {alpha_err:.3f}")
    except Exception as e:
        print(f"Error fitting Rg vs N: {e}")

def main():
    # Configure these parameters based on your simulation setup
    base_dir = "analysis"  # Changed from "dumps" to "analysis"
    output_dir = "results"
    n_values = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]  # Example N values
    seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]  # Example seed values
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process MSD data
    print("Processing MSD data...")
    msd_results = process_msd_data(base_dir, n_values, seeds)
    
    # Check if any MSD data was found
    if not msd_results:
        print("No MSD data found! Check file paths and directory structure.")
        return
    
    # Process Rg data
    print("Processing radius of gyration data...")
    rg_results = process_rg_data(base_dir, n_values, seeds)
    
    # Check if any Rg data was found
    if not rg_results:
        print("No radius of gyration data found! Check file paths and directory structure.")
        return
    
    # Calculate diffusion coefficients
    print("Calculating diffusion coefficients...")
    diff_coeff = calculate_diffusion_coefficients(msd_results)
    
    # Create plots
    print("Creating plots...")
    if msd_results:
        plot_msd_vs_time(msd_results, output_dir)
    if diff_coeff:
        plot_diffusion_vs_n(diff_coeff, output_dir)
    if rg_results:
        plot_rg_vs_n(rg_results, output_dir)
    
    print(f"Analysis complete. Results saved to {output_dir}")

if __name__ == "__main__":
    main()