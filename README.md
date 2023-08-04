# PARROT (**P**rotein allocation **A**djustment fo**R** alte**R**native envir**O**nmen**T**s)
Constraint-based prediction of protein abundance for alternative conditions

## Dependencies
* MATLAB (tested with 2017a, 2022a and 2023a)
* [COBRA toolbox 3.0](https://github.com/opencobra/cobratoolbox)
* [RAVEN toolbox 2.0](https://github.com/SysBioChalmers/RAVEN)
* [GECKO toolbox 2.0](https://github.com/SysBioChalmers/GECKO)
* [Gurobi solver](https://www.gurobi.com/) (tested with version 9.1.1)

## How to run PARROT for your model
1. Create a pcGEM using the [GECKO toolbox](https://github.com/SysBioChalmers/GECKO)
2. Integrate protein abundance values for the reference condition in the created pcGEM 
3. Define the constraints and simulation parameters for the alternative condition
4. Run PARROT to obtain the enzyme usage distribution for your target alternative condition

| Input | Explanation |
| :---:         | --- |
| _modelREF_     | Protein-constrained model with integrated protein measurements to be used as a reference condition. |
| _modelSTR_     | Protein-constrained model, batch version. |
| _MinStr_   | Minimization strategy, 'Euclidean' or 'Manhattan' |
| _lambda_   | Weighting factor for metabolic fluxes |
| _etotREF_   | Total enzyme usage for the reference condition |
| _etotALT_   | Total enzyme usage for the alternative condition |

## Usage:
[solutionALT, solStatus] = PARROT(modelREF, modelALT, MinStr, lambda, etotREF, etotALT)

Refer to 'Example1.m' and 'Example2.m' for additional information on how to run PARROT.

