# UF6 Hydrolysis Kinetics Estimator

## Motivation
Uranium hexafluoride, UF6, is a critical compound in the nuclear fuel cycle due to its use in uranium enrichment. 
Despite its high importance, UF6 presents significant safety challenges. 
When UF6 makes contact with moisture, it hydrolyzes to produce radioactive uranyl fluoride particulates, and acidic hydrogen fluoride gas, both of which are hazardous to human health and the environment. 
Accurate prediction of UF6 behavior requires reliable knowledge regarding its reaction kinetics; however, the hydrolysis process remains difficult to decipher. 
The hydrolysis process happens rapidly making it difficult to understand the entirety of the reaction.
Improving kinetics understanding is therefore essential for national security, computational modeling, and development of effective mitigation strategies. 

Building from previous research, a new model where UF6 reacts with a water dimer rather than seperate molecules is the basis of new research. 
This repository attempts to advance the current study of UF6 Hydrolysis for researchers by providing scientific software that implements the new kinetic model for future scientific research. 

## Files

- `effectivek_estimator.py`: This software is the main focus of the repository. This script takes in methods and parameters for calculating the effective rate constant from user inputs and returns the effective rate constant.
- `hydrolysis.py`: This script calculates the rate of reaction for the kinetic model and compares it to values from experimental rates and rates from previous models. 
- `k_optimize.py`: This script utilizes the Pyomo package to optimize the parameter space of the Hertz-Knudsen pressure decay model.
- `k_optimize_two.py`: This script implements a manual grid search to optimize the parameter space of the Hertz-Knudsen pressure decay model.
- `rateconst.py`: This script calculates the rate constants of desoprtion and adsorpton given 

## Logging:
Please note that the current logging level is set to INFO. If you wish to change this, open `effectivek_estimator.py` with a text or code editor and replace the INFO in line 5 to whichever level you want to run. (DEBUG, INFO, WARNING, ERROR, CRITICAL) 

## Usage
