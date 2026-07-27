import numpy as np
import matplotlib.pyplot as plt
import math

'''
This Python file conducts a manual grid search/optimization of the Hertz-Knudsen approximization to the rate of reaction of UF6

This problem contains various variable parameters that are within bounds. These bounds act as physical constraints however, more physical constraints can and should be added for a more accurate model.
'''


def build_data():
    # the data dictionary contains the partial pressures of H2O and UF6 as well as the corresponding observed rate constant for each Richards experiment

    data = {
        1: {'PH2O': 60.0, 'PUF6': 8.2,  'kobs': 0.04719},
        2: {'PH2O': 60.0, 'PUF6': 16.2, 'kobs': 0.03432},
        3: {'PH2O': 60.0, 'PUF6': 30.2, 'kobs': 0.02347},
        4: {'PH2O': 80.0, 'PUF6': 9.1,  'kobs': 0.07714},
        5: {'PH2O': 80.0, 'PUF6': 18.6, 'kobs': 0.05473},
        6: {'PH2O': 80.0, 'PUF6': 28.3, 'kobs': 0.04636}
    }
    return data

def term_two(T = 298.0, m = 5.845e-25, kB = 1.380649e-23):
    # this function returns the second term of the k_eff approx for readable computation

    return math.sqrt(kB * T / (2.0 * math.pi * m))

def water_available_fraction(PH2O_mTorr, T = 298.0):
    '''
    This function finds the fraction of available water for the water loading constraint

    Given some partial pressure of water introduced in the experiment, this function attempts to utilize it as a number density constraint, be returning the available water fraction for a given water vapor pressure and temperature
    fraction
    '''

    P = PH2O_mTorr * 0.133322   # Pa
    Rgas = 8.314462618          # J/mol/K
    Mw = 0.01801528             # kg/mol
    rho = 1000.0                # kg/m^3

    phi_avail = P * Mw / (Rgas * T * rho)
    return phi_avail

def generate_tradeoff(kobs, PH2O, R_bounds, K_bounds):
    
    """
    This function generates the feasible parameter combinations for a single experiment and adds the result to a dictionary. The number of points used is manually adjustable in the function 

    Additional constraints can be added to steer away from degeneracy.
    """

    alpha = 1.0
    T = 298.0
    m = 5.845e-25
    kB = 1.380649e-23
    nR = 150
    nK = 150

    # calculate second term of approximation
    term2 = term_two(T = T, m = m, kB = kB)

    # create a grid of plausible Radius and K values based on the bounds and number of points in the parameter space. 
    R_grid = np.logspace(np.log10(R_bounds[0]), np.log10(R_bounds[1]), nR)
    K_grid = np.logspace(np.log10(K_bounds[0]), np.log10(K_bounds[1]), nK)

    # creating a "mesh" indexed by radius and K 
    R_mesh, K_mesh = np.meshgrid(R_grid, K_grid, indexing='ij')

    # create all plausible ranges of number density based on the Radius and K grid that exactly satisfies the k_eff equation
    nd_mesh = kobs / (4.0 * np.pi * R_mesh**2 * term2 * alpha * K_mesh)

    # ----------- need more constraints for the number density or anything in general ---------
    # liquid loading / volume fraction 
    # use the available PH2O from the data... 
    # adding constraint
    phi_max = water_available_fraction(PH2O)
    loading_mesh = nd_mesh * (4.0 / 3.0) * np.pi * R_mesh**3
    load_constraint = loading_mesh <= phi_max

    filter = load_constraint

    # add filter and return results
    results = {
    'filter': filter,
    'R_mesh': R_mesh,
    'K_mesh': K_mesh,
    'nd_mesh': nd_mesh,
    'loading_mesh': loading_mesh,
    'R_valid': R_mesh[filter],
    'K_valid': K_mesh[filter],
    'nd_valid': nd_mesh[filter],
    'loading_valid': loading_mesh[filter]
    }

    return results


def plot_tradeoffs(exp_id, exp_data, results):

    '''
    This function recieves a experiment id and its data/results from the tradeoff function for plotting.

    Each plot contains the combination of droplet number density, droplet radius, and the K fraction that fit the k_eff equation and constraints
    '''

    K_valid = results['K_valid']
    R_valid = results['R_valid']
    nd_valid = results['nd_valid']

    if len(R_valid) == 0:
        print(f"No \feasible points for experiment {exp_id}")
        return

    plt.style.use('seaborn-v0_8-whitegrid')

    # Generate subplots for clean visualization

    fig, axs = plt.subplots(1, 3, figsize=(20, 6))

    fig.suptitle(
        f"Experiment {exp_id}: "
        f"PH2O={exp_data['PH2O']} mTorr, "
        f"PUF6={exp_data['PUF6']} mTorr, "
        f"kobs={exp_data['kobs']:.5f} s$^{{-1}}$",
        fontsize=12
    )

    # nd vs R
    sc1 = axs[0].scatter(
        R_valid, nd_valid, c=K_valid, cmap='viridis', s=18, alpha=0.8
    )
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('Droplet radius, R (m)')
    axs[0].set_ylabel(r'Droplet density, $N_d/V$ (m$^{-3}$)')
    axs[0].set_title(r'Tradeoff: $N_d/V$ vs R')
    fig.colorbar(sc1, ax=axs[0], label='K')

    # nd vs K
    sc2 = axs[1].scatter(
        K_valid, nd_valid, c=R_valid, cmap='plasma', s=18, alpha=0.8
    )
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Reactive probability, K')
    axs[1].set_ylabel(r'Droplet density, $N_d/V$ (m$^{-3}$)')
    axs[1].set_title(r'Tradeoff: $N_d/V$ vs K')
    fig.colorbar(sc2, ax=axs[1], label='R (m)')

    # K vs R
    sc3 = axs[2].scatter(
        R_valid, K_valid, c=nd_valid, cmap='cividis', s=18, alpha=0.8
    )
    axs[2].set_xscale('log')
    axs[2].set_yscale('log')
    axs[2].set_xlabel('Droplet radius, R (m)')
    axs[2].set_ylabel('Reactive probability, K')
    axs[2].set_title('Tradeoff: K vs R')
    fig.colorbar(sc3, ax=axs[2], label='Droplet density')


    plt.tight_layout()
    plt.show()


def print_ranges(exp_id, results):
    '''
    This function takes the experiment id and its results from the tradeoff function, and prints the range of values into the terminal
    '''

    if len(results['R_valid']) == 0:
        print(f"Experiment {exp_id}: no feasible solutions")
        return

    print(f"\nExperiment {exp_id}: feasible parameter ranges")
    print(f"R:  {results['R_valid'].min():.3e} to {results['R_valid'].max():.3e} m")
    print(f"K:  {results['K_valid'].min():.3e} to {results['K_valid'].max():.3e}")
    print(f"Nd/V: {results['nd_valid'].min():.3e} to {results['nd_valid'].max():.3e} m^-3")


def plot_fit_all(global_results, tolerance):
    '''
    This function attempts to combine the optimization results of all six experiments based on some given tolerance

    Radius and K values are set the same for each experiment however the number density varies and none of them are the same 

    Using the tolerance we can find some maximum multiple of the miimum density value that would be considered for plotting
    '''

    R = global_results[0]['R_mesh']
    K = global_results[0]['K_mesh']

    # Stack number density arrays
    nd_stack = np.stack(
        [exp['nd_mesh'] for exp in global_results],
        axis=0
    )

    nd_mean = np.mean(nd_stack, axis=0)
    nd_min = np.min(nd_stack, axis=0)
    nd_max = np.max(nd_stack, axis=0)
    nd_ratio = nd_max / nd_min

    # Require number densities across experiments to be close not exact within some ratio tolerance 
    close_filter = nd_ratio <= tolerance

    fit_filter = close_filter

    if not np.any(fit_filter):
        print("No points found. Change tolerance.")
        return

    plt.figure(figsize=(8, 6))

    # begin scatter plot with the average number density across experiments as the Y and radius as the x, highlighting the K values through the plot
    sc = plt.scatter(
        R[fit_filter],
        nd_mean[fit_filter],
        c=K[fit_filter],
        cmap='viridis',
        s=30,
        alpha=0.8
    )

    plt.colorbar(sc, label='K')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Droplet radius, R (m)')
    plt.ylabel(r'Mean droplet density, $N_d/V$')
    plt.title(f'Experimental fit region: Average Number Droplet Density vs Radius, Tolerance={tolerance}'
    )

    plt.tight_layout()
    plt.show()

    # print results
    print("Fit region:")
    print(f"Mean Nd/V: {nd_mean[fit_filter].min():.3e} to {nd_mean[fit_filter].max():.3e} m^-3")
    print(f"Nd ratio: {nd_ratio[fit_filter].min():.3f} to {nd_ratio[fit_filter].max():.3f}")

def main():
    data = build_data()

    # plausible ranges for the model
    R_bounds = (1e-8, 1e-6)
    K_bounds = (0.45, 0.55)

    # plot for all experiments 

    global_results = []

    for exp_id, exp_data in data.items():

        pressure = data[exp_id]['PH2O']

        results = generate_tradeoff(
            kobs=exp_data['kobs'],
            PH2O = pressure,
            R_bounds=R_bounds,
            K_bounds=K_bounds
        )
        global_results.append(results)
        print_ranges(exp_id, results)
        plot_tradeoffs(exp_id, exp_data, results)

    # global plot 

    plot_fit_all(global_results, tolerance=3.5)

if __name__ == "__main__":
    main()