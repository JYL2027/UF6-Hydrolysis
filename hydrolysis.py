import numpy as np
import matplotlib.pyplot as plt

def keq(MESS: str) -> list:
    '''
    This code takes in a input file (MESS.out) and calculates the keq rate constant at temperature T at the high pressure limit

    Parameters: A input file, MESS.out

    Output: A list of keq rate constants
    '''
    # capture/escape rates vs T @ High Pressure
    # Store arrays for T and rates
   
    Start_Process = False
    Temp = []
    ls1_rate = []
    rs1_rate = []
    b1_rate = []
    keq = []

    try:
        with open(MESS, 'r') as file:
            for line in file:
                if line.strip() == "Capture/Escape Rate Coefficients:":
                    Start_Process = True
                    continue

                if Start_Process:
                    if line.strip().startswith("______________________________________________________________________________________"):
                        break
                    parts = line.split()

                    if len(parts) < 1 or parts[0] == "T(K)":
                            continue
                    
                    Temp.append(float(parts[0]))
                    ls1_rate.append(float(parts[1]))
                    rs1_rate.append(float(parts[2]))
                    b1_rate.append(float(parts[3]))
    except Exception as e:
        print(f"Error reading {MESS}: {e}")
        return None

    # Calculate keq rate constant
    # Based on file readings it seems that this value is equal to b_1/ls_1
    # cm^3/molecule

    try: 
        for i in range(len(ls1_rate)):
            keq.append(b1_rate[i]/ls1_rate[i])
    except Exception as e:
        print(f"Error in calculations: {e}")
        return None
    
    return keq

def k_2(MESS: str) -> list:
    '''
    This function takes the MESS output file and finds the forward reaction rate k2 illustrated in the paper.

    Parameters: A input file, MESS.out

    Output: A list of keq rate constants
    '''

    Start_Process = False
    k2 = []

    try:
        with open(MESS, 'r') as file:
            for line in file:
                if line.strip().startswith("Reactant = ls_"):
                    Start_Process = True
                    continue

                if Start_Process:
                    if line.strip().startswith("Reactant = rs_"):
                        break
                    parts = line.split()

                    if len(parts) < 1 or parts[0] == "T(K)":
                            continue
                    
                    k2.append(float(parts[2]))
    except Exception as e:
        print(f"Error reading {MESS}: {e}")
        return None

    return k2

def waterkeq(MESS: str) -> list:
   
    Start_Process = False
    Temp = []
    ls1_rate = []
    rs1_rate = []
    b1_rate = []
    keq = []

    try:
        with open(MESS, 'r') as file:
            for line in file:
                if line.strip() == "Capture/Escape Rate Coefficients:":
                    Start_Process = True
                    continue

                if Start_Process:
                    if line.strip().startswith("______________________________________________________________________________________"):
                        break
                    parts = line.split()

                    if len(parts) < 1 or parts[0] == "T(K)":
                            continue
                    
                    Temp.append(float(parts[0]))
                    ls1_rate.append(float(parts[1]))
                    b1_rate.append(float(parts[2]))
    except Exception as e:
        print(f"Error reading {MESS}: {e}")
        return None

    try: 
        for i in range(len(ls1_rate)):
            keq.append(b1_rate[i]/ls1_rate[i])
    except Exception as e:
        print(f"Error in calculations: {e}")
        return None
    
    return keq

def calc_droplet_number_density(k_eff, radius):

    T=295.0
    m=5.845e-25
    alpha=1.0
    reactive_fraction=0.5

    k_b = 1.380649e-23

    thermal_factor = np.sqrt(k_b * T / (2.0 * np.pi * m))

    area = 4.0 * np.pi * radius**2

    droplet_nd = k_eff / (area * thermal_factor * alpha * reactive_fraction)

    return droplet_nd

def main():
    new = ["CC2_UF6H2OH2Oa.out", "CC2_UF6H2OH2Ob.out"]
    dimers = ["CC2_UF6H2OH2Oa.out", "CC2_UF6H2OH2Ob.out"]

    # Calculating the effective rate and implementing the calculations of the rate of reaction
    # Concentrations from Richards in order of uf6 concentration
    CH2O = [60, 80, 60, 80, 80, 60]
    CUF6 = [8.2, 9.1, 16.2, 18.6, 28.3, 30.2]
    # Reaction rates of the concentrations from Richards
    Exp_rate = [0.387, 0.702, 0.556, 1.018, 1.312, 0.709]
    Exp_error = [0.008, 0.1, 0.071, 0.066, 0.112, 0.321]

    # Convert to numpy arrays
    CH2O = np.array(CH2O)
    CUF6 = np.array(CUF6)
    Exp_rate = np.array(Exp_rate)
    Exp_error = np.array(Exp_error)

    # Constants
    R = 0.082057          # L atm mol^-1 K^-1
    T = 295.0             # K
    N_A = 6.02214076e23   # molecules/mol
    # Convert mTorr -> Torr
    CH2O_torr = CH2O / 1000.0
    CUF6_torr = CUF6 / 1000.0
    # Convert Torr -> atm
    CH2O_atm = CH2O_torr / 760.0
    CUF6_atm = CUF6_torr / 760.0
    # Convert atm -> mol/L
    CH2O_M = CH2O_atm / (R * T)
    CUF6_M = CUF6_atm / (R * T)
    # Convert mol/L -> molecules/cm^3
    CH2O_numdens = CH2O_M * N_A / 1000.0
    CUF6_numdens = CUF6_M * N_A / 1000.0

    Exp_rate = Exp_rate / 1000 / 760.0 / (R * T)
    Exp_error = Exp_error / 1000 / 760.0 / (R * T)
    # -------------------------------------------------------

    # Create dictionaries for the two rates 
    adsorption = {
        1e-8:  9.991e9,
        1e-7: 9.991e8,
        1e-6: 9.991e7
    }

    desorption = {
        (30, 4.397e12): 2.143e+07, 
        (30, 1.112e13): 5.419e+07, 
        (30, 1.352e13): 6.589e+07, 
        (40, 4.397e12): 3.633e5, 
        (40, 1.112e13): 9.187e5, 
        (40, 1.352e13): 1.117e6, 
        (45, 4.397e12): 4.730e4, 
        (45, 1.112e13): 1.196e5, 
        (45, 1.352e13): 1.454e5, 
        (55, 4.397e12): 8.019e+02, 
        (55, 1.112e13): 2.028e+03, 
        (55, 1.352e13): 2.466e+03
    }

    # Start figures
    for radius_val, k_ads in adsorption.items():
        for (E_des, nu), k_des in desorption.items():

            plt.figure(figsize=(8, 6))

            plt.xlabel("Concentration of UF6 (M)")
            plt.ylabel("Rate of Reaction (M/s)")
            plt.title(
                f"Experimental vs Computational Rates\n"
                f"radius = {radius_val}, "
                f"$k_{{ads}}$ = {k_ads:.3e} s$^{{-1}}$, "
                f"$E_{{des}}$ = {E_des} kJ/mol, "
                f"$\\nu$ = {nu:.3e} s$^{{-1}}$, "
                f"$k_{{des}}$ = {k_des:.3e} s$^{{-1}}$"
            )

            plt.errorbar(
                CUF6_M, Exp_rate, yerr=Exp_error,
                fmt='-x', capsize=5, ecolor='black',
                label='Experiment'
            )

            # code the original outputs from the paper
            for dimerfile in dimers:
                k_dime = waterkeq("CC_H2O.out")
                k_constant_eq = keq(dimerfile)
                k2 = k_2(dimerfile)

                Comp_Rates = []
                for idx in range(len(CH2O_numdens)):
                    rate = (
                        k_dime[0]
                        * k_constant_eq[0]
                        * k2[0]
                        * CUF6_numdens[idx]
                        * CH2O_numdens[idx]
                        * CH2O_numdens[idx]
                    )
                    Comp_Rates.append(rate)

                Comp_Rates = np.array(Comp_Rates)
                Comp_Rates = Comp_Rates * 1000 / N_A

                plt.plot(CUF6_M, Comp_Rates, '-o', label=f'{dimerfile}')


            # use the new hypothesis for graphing with the preliminary model, including number density of the dropelts at each point on the plot
            for filer in new:
                nd_val = []
                k_constant_eq = keq(filer)
                k2 = k_2(filer)
                
                Comp_Rates = []
                for idx in range(len(CH2O_numdens)):

                    # the k_eff is k_ads / k_des * k_constant_eq[0] * k2[0] * CH2O_numdens[idx] corresponding to the model of the Hertz-Knudsen model dP/dt = k_eff[PUF6]

                    nd = calc_droplet_number_density((k_ads / k_des * k_constant_eq[0] * k2[0] * CH2O_numdens[idx]), radius = radius_val)

                    nd_val.append(nd)
                    rate = (
                        k_ads / k_des * k_constant_eq[0] * k2[0] * CUF6_numdens[idx] * CH2O_numdens[idx]
                    )
                    Comp_Rates.append(rate)

                Comp_Rates = np.array(Comp_Rates)
                Comp_Rates = Comp_Rates * 1000 / N_A

                plt.plot(
                    CUF6_M,
                    Comp_Rates,
                    '-o',
                    label=f'{filer}.update'
                )

                for i, label in enumerate(nd_val):
                    plt.annotate(
                        f"{label:.1e}",                    
                        xy=(CUF6_M[i], Comp_Rates[i]),            
                        xytext=(8, 8),              
                        textcoords='offset points'  
                    )
                    
            plt.yscale('log')
            plt.grid(True)
            plt.legend(fontsize="small")
            plt.tight_layout()
            plt.show()

if __name__ == "__main__":
    main()