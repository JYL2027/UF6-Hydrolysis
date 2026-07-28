# UF6 Hydrolysis Kinetics Estimator

## Objective:
Uranium hexafluoride, UF6, is a critical compound in the nuclear fuel cycle due to its use in uranium enrichment. 
Despite its high importance, UF6 presents significant safety challenges. 
When UF6 makes contact with moisture, it hydrolyzes to produce radioactive uranyl fluoride particulates, and acidic hydrogen fluoride gas, both of which are hazardous to human health and the environment. 
Accurate prediction of UF6 behavior requires reliable knowledge regarding its reaction kinetics; however, the hydrolysis process remains difficult to decipher. 
The hydrolysis process happens rapidly making it difficult to understand the entirety of the reaction.
Improving kinetics understanding is therefore essential for national security, computational modeling, and development of effective mitigation strategies. 

Building from previous research, a new model where UF6 reacts with a water dimer rather than seperate molecules is the basis of new research. 
This repository attempts to advance the current study of UF6 Hydrolysis for researchers by providing scientific software that implements the new kinetic model for future scientific research. 

## Files:

- `effectivek_estimator.py`: This software is the main focus of the repository. This script takes in methods and parameters for calculating the effective rate constant from user inputs and returns the effective rate constant.
- `hydrolysis.py`: This script calculates the rate of reaction for the kinetic model and compares it to values from experimental rates and rates from previous models. 
- `k_optimize.py`: This script utilizes the Pyomo package to optimize the parameter space of the Hertz-Knudsen pressure decay model.
- `k_optimize_two.py`: This script implements a manual grid search to optimize the parameter space of the Hertz-Knudsen pressure decay model.
- `rateconst.py`: This script calculates the rate constants of desoprtion and adsorpton given a variety of parmeters.

## Logging:
Please note that the current logging level is set to INFO. If you wish to change this, open `effectivek_estimator.py` with a text or code editor and replace the INFO in line 5 to whichever level you want to run. (DEBUG, INFO, WARNING, ERROR, CRITICAL) 

## Usage:

To use the `effectivek_estimator.py` software please first clone this repository. From there you can run the script in whichever Python interpretor you wish. 

# Calculation methods:
The software has three methods for calulcating the effective rate constant. 
- `MESS_HK`: This method uses a MESS based approach that utilizes the Hertz-Knudsen to approximate k_ads. 
This method returns two effective rate constants. One in the form `rate = keff[H2O][UF6]` and one in the form `rate = keff[UF6]`. 
(`keff = keq * k_2 * kads/kdes or keff = keq * k_2 * kads/kdes * [H2O]`)
- `HK_P`: This method produces the Hertz-Knudsen pressure decay effective rate constant of the form `dP/dt = keff [PUF6]`. 
- `BEMM`: This method produces a pressure decay effective rate constant using the Bolzmann Equation Moment Method model as a estimate to adsorpton. 
- `LinearBEMM`: This method produces a pressure decay effective rate constant using the linearized Bolzmann Equation Moment Method model.
# Parameters:
Each method available has a unique set of input parameters that must be met for the solver to calculate the effective rate constant. Please enter the parameters seperated by commas in the following order for each method.
- `MESS_HK`: <br> Equilibrium rate constant for dissociation from MESS (keq), <br> Forward reaction rate constant from MESS (k2), <br> Radius of a water micro droplet (m), <br> Accommodation Coefficient, <br> Temperature of reaction (K), <br> Attempt frequency (1/s), <br> Desorption Energy (kJ/mol), <br> Concentration of water dimers (molecules/cm^3)
- `HK_P`:  <br> Radius of a water micro droplet (m),<br> droplet number density (droplets/m^3), <br> Energy for reaction (kJ/mol),<br> Energy for deesorption (kJ/mol),<br>Temperature of reaction (K), <br> The attempt frequency (1/s), <br>Accommodation Coefficient
- `BEMM`: <br> Pressure of the UF6 gas/vapor (Pa), <br> Temperature of the surface of the water-micro-droplet (K), <br> Temperature of the UF6 vapor/gase (K), <br> Velocity of the UF6 gas molecules (m/s), <br> Radius of the water-microdroplet (m), <br> Droplet number density (droplets/m^3), <br> The probability of reaction before desoprtion
- `LinearBEMM`: <br> Alpha the accomodation coefficient, <br> The temperature of the water micro droplet surface (K), <br> Radius of the micro droplet (m),
        <br> Temperature of the UF6 gas/vapor (K), 
        <br> Pressure of the Uf6 gas/vapor (Pa),
        <br> Effective UF6 partial pressure at the water droplet surface (Pa),
        <br> The droplet number density of the water micro droplet (droplets/m^3),
        <br> The Probabilty of reaction before desorbing

# Example Execution:

Code: 
```
Available methods:
MESS_HK
HK_P
BEMM
LinearBEMM
What method would you like to use?
```

User:
```
HK_P
```

Code:
```
Please input all the data necessary to compute the effective rate constant. Please enter comma-separated numeric values.
```

User:
```
1e-6, 1e12, 25, 40, 298, 1e13, 0.1
```

Code:
```
The effective rate constant or constants with method HK_P is/are [41.96227663622824] 1/s
```

