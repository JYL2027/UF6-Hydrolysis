import pyomo.environ as pyo
import math
import matplotlib.pyplot as plt
import numpy as np

def build_data():

    # load data as a dictionary. 
    # rate constants of each experiments are kobs, while partial pressure of the reactants are included
    data = {
        1: {'PH2O': 60.0, 'PUF6': 8.2,  'kobs': 0.04719},
        2: {'PH2O': 60.0, 'PUF6': 16.2, 'kobs': 0.03432},
        3: {'PH2O': 60.0, 'PUF6': 30.2, 'kobs': 0.02347},
        4: {'PH2O': 80.0, 'PUF6': 9.1,  'kobs': 0.07714},
        5: {'PH2O': 80.0, 'PUF6': 18.6, 'kobs': 0.05473},
        6: {'PH2O': 80.0, 'PUF6': 28.3, 'kobs': 0.04636}
    }

    temperature = 298
    boltz = 1.380649e-23
    mass = 5.845e-25
    gas_const = 0.008314

    return data, temperature, boltz, mass, gas_const


def build_model(data, temperature, boltz, mass, gas_const, target_exp):

    m = pyo.ConcreteModel()

    # index everyting be experiment number
    m.exp = pyo.Set(initialize=data.keys())

    m.PH2O = pyo.Param(m.exp, initialize={i: data[i]['PH2O'] for i in data})
    m.PUF6 = pyo.Param(m.exp, initialize={i: data[i]['PUF6'] for i in data})
    m.kobs = pyo.Param(m.exp, initialize={i: data[i]['kobs'] for i in data})

    # establishing constant parameters
    m.Temp = pyo.Param(initialize=temperature)
    m.Boltz = pyo.Param(initialize=boltz)
    m.Mass = pyo.Param(initialize=mass)
    m.Gas_Const = pyo.Param(initialize=gas_const)

    # Establishing variables for the model
    # keep nu and alpha constant for now as more variables is worse for the model
    m.nu = pyo.Param(initialize=1.2e13)
    m.alpha = pyo.Param(initialize=1.0)

    # we want the radius, number density, and K ratio to be parameters adjustable to the model
    m.R = pyo.Var(bounds=(1e-8, 1e-6), initialize=1e-7, within=pyo.PositiveReals)
    m.nd = pyo.Var(initialize=1e10, within=pyo.PositiveReals)
    m.K = pyo.Var(bounds=(1e-4, 1.0), initialize=0.5, within=pyo.PositiveReals)

    # ------------Adding constraints--------------------
    # liquid loading constraint suggested by AI
    # Plan to fix K so need a relationship fo nd and R

    phi_max = 1e-6
    m.liquid = pyo.Constraint(
        expr=m.nd * (4.0 / 3.0) * math.pi * m.R**3 <= phi_max
    )

    # -----------Creating Expressions and objectives--------
    m.thermal = pyo.Expression(
        expr=(m.Boltz * m.Temp / (2.0 * math.pi * m.Mass))**0.5
    )

    m.kpred = pyo.Expression(
        expr=4.0 * math.pi * m.R**2 * m.nd * m.thermal * m.alpha * m.K
    )

    # Objective for one experiment only
    # minimize the squared absolute eror
    m.obj = pyo.Objective(
        expr=(m.kpred - m.kobs[target_exp])**2,
        sense=pyo.minimize
    )

    return m


def solve(model, Kval, target_exp, solver_name='ipopt'):
    """
    This function fixes K values, and solve the model for nd and R for one target experiment.
    """

    model.K.fix(Kval)

    solver = pyo.SolverFactory(solver_name)
    results = solver.solve(model, tee=False)

    out = {
        'exp': target_exp,
        'R': pyo.value(model.R),
        'K': Kval,
        'nd': pyo.value(model.nd),
        'obj': pyo.value(model.obj),
        'kpred': pyo.value(model.kpred),
        'kobs': pyo.value(model.kobs[target_exp]),
    }

    out['rel_err'] = abs(out['kpred'] - out['kobs']) / out['kobs']

    model.K.unfix()

    return out


def plot_solutions(solutions, exp_id, data):
    """
    Create tradeoff plots for one experiment, these plots are scatter plots and create various comparisons. Each dot is colored based on their relative error. 
    """

    if len(solutions) == 0:
        print(f"No solutions to plot for experiment {exp_id}.")
        return

    R_vals = [sol['R'] for sol in solutions]
    nd_vals = [sol['nd'] for sol in solutions]
    K_vals = [sol['K'] for sol in solutions]
    obj_vals = [sol['obj'] for sol in solutions]
    err_vals = [sol['rel_err'] for sol in solutions]

    PH2O = data[exp_id]['PH2O']
    PUF6 = data[exp_id]['PUF6']
    kobs = data[exp_id]['kobs']

    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        f'Experiment {exp_id}: PH2O={PH2O} mTorr, PUF6={PUF6} mTorr, kobs={kobs:.5f} s^-1',
        fontsize=13
    )

    # nd vs R
    sc1 = axs[0, 0].scatter(R_vals, nd_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[0, 0].set_xscale('log')
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_xlabel('R (m)')
    axs[0, 0].set_ylabel('nd (1/m^3)')
    axs[0, 0].set_title('nd vs R')
    axs[0, 0].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc1, ax=axs[0, 0], label='Relative Error')

    # nd vs K
    sc2 = axs[0, 1].scatter(K_vals, nd_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_yscale('log')
    axs[0, 1].set_xlabel('K')
    axs[0, 1].set_ylabel('nd (1/m^3)')
    axs[0, 1].set_title('nd vs K')
    axs[0, 1].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc2, ax=axs[0, 1], label='Relative Error')

    # K vs R
    sc3 = axs[0, 2].scatter(R_vals, K_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[0, 2].set_xscale('log')
    axs[0, 2].set_yscale('log')
    axs[0, 2].set_xlabel('R (m)')
    axs[0, 2].set_ylabel('K')
    axs[0, 2].set_title('K vs R')
    axs[0, 2].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc3, ax=axs[0, 2], label='Relative Error')

    # objective vs nd
    sc4 = axs[1, 0].scatter(nd_vals, obj_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[1, 0].set_xscale('log')
    axs[1, 0].set_xlabel('nd (1/m^3)')
    axs[1, 0].set_ylabel('Objective')
    axs[1, 0].set_title('Objective vs nd')
    axs[1, 0].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc4, ax=axs[1, 0], label='Relative Error')

    # objective vs R
    sc5 = axs[1, 1].scatter(R_vals, obj_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_xlabel('R (m)')
    axs[1, 1].set_ylabel('Objective')
    axs[1, 1].set_title('Objective vs R')
    axs[1, 1].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc5, ax=axs[1, 1], label='Relative Error')

    # objective vs K
    sc6 = axs[1, 2].scatter(K_vals, obj_vals, c=err_vals, cmap='viridis', alpha=0.85)
    axs[1, 2].set_xscale('log')
    axs[1, 2].set_xlabel('K')
    axs[1, 2].set_ylabel('Objective')
    axs[1, 2].set_title('Objective vs K')
    axs[1, 2].grid(True, which='both', ls='--', alpha=0.4)
    fig.colorbar(sc6, ax=axs[1, 2], label='Relative Error')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def main():
    data, temperature, boltz, mass, gas_const = build_data()

    # R_grid = np.logspace(-8, -6, 15)
    K_grid = np.logspace(-4, 0, 50)

    all_filtered = []
 
    # Solve separately for each experiment
    for exp_id in data.keys():
        print(f"\nRunning experiment {exp_id}...")

        solutions = []
        m = build_model(data, temperature, boltz, mass, gas_const, exp_id)
        for Kval in K_grid:
                result = solve(m,Kval, exp_id)
                solutions.append(result)

        # filter "family" based on some rel error threshhold.
        filtered = [sol for sol in solutions if sol['rel_err'] < 0.1]

        print(f"Experiment {exp_id}: total solutions = {len(solutions)}")
        print(f"Experiment {exp_id}: filtered solutions = {len(filtered)}")

        if len(filtered) > 0:
            print("R min/max:", min(sol['R'] for sol in filtered), max(sol['R'] for sol in filtered))
            print("K min/max:", min(sol['K'] for sol in filtered), max(sol['K'] for sol in filtered))
            print("nd min/max:", min(sol['nd'] for sol in filtered), max(sol['nd'] for sol in filtered))
            plot_solutions(filtered, exp_id, data)
            all_filtered.extend(filtered)
        else:
            print(f"No filtered family for experiment {exp_id}; plotting all solutions.")
            plot_solutions(solutions, exp_id, data)
            all_filtered.extend(solutions)

 # This still does not fully identify R and nd .


if __name__ == "__main__":
    main()