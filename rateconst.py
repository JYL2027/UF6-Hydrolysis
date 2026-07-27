import numpy as np

# This code calculates the adsorption and desorption constants for a variety of inputs. 
Boltz = 1.380649e-23
T = 295 # K
p = 1.999 #Pa
m = 5.845e-25 #kg
alpha = 1
radius = [1e-8, 1e-7, 1e-6] #m
energy = [30, 40, 45, 55] # kJ/mol
frequency = [4.397e12, 1.112e13, 1.352e13] # 1/s
R = 0.008314 

def adsorption_rate(radius):
    k = (p * alpha) / (np.sqrt(2*np.pi*m*Boltz*T))
    k_ads = k * 4 * np.pi * radius**2

    # Convert to a rate constant

    k_rate_const = k_ads * (1/(p/(Boltz*T))) * (1/((4/3)*np.pi*((radius)**3)))

    return k_rate_const

# Desorption calculations 

def desorption_rate(energy, frequency):

    return frequency * np.exp(-energy / (R * T))

def main():
    # Find adsorption rates
    for r in radius: 
        rate = adsorption_rate(r)
        print(f"The adsorption rate constant for alpha = 1, radius = {r} is {rate:.3e} 1/s")

    # Find desorption rates
    for i in energy:
        for j in frequency:
            rate = desorption_rate(i, j)
            print(f"The desorption rate constant for a Desorption Energy of {i} kJ/mol and attempt frequency of {j:.3e} 1/s is {rate:.3e} 1/s")

if __name__ == "__main__":
    main()