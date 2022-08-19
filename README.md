# PARROT (**P**rotein allocation **A**djustment fo**R** st**R**ess c**O**ndi**T**ions)
Constraint-based prediction of protein abundance for suboptimal conditions

## Dependencies
* MATLAB (tested with 2017a and 2022a)
* [COBRA toolbox 3.0](https://github.com/opencobra/cobratoolbox)
* [RAVEN toolbox 2.0](https://github.com/SysBioChalmers/RAVEN)
* [GECKO toolbox 2.0](https://github.com/SysBioChalmers/GECKO)
* [Gurobi solver](https://www.gurobi.com/) (tested with version 9.1.1)

## How to run PARROT for your model
1. Create an ecModel using the [GECKO toolbox](https://github.com/SysBioChalmers/GECKO)
2. Integrate protein abundance values for the reference condition in the created ecModel 
3. Adjust in the PARROT code the constraints that correspond to your target suboptimal condition
4. Adjust the file names, file paths and relevant parameters 
5. Run PARROT to obtain the enzyme usage distribution for your target suboptimal condition
the weighted 

| Input | Explanation |
| :---:         | --- |
| _modelREF_     | pcGEM integrated with experimental proteomics for reference condition |
| _modelSTR_   | pcGEM without integrated proteomics data |
| _modelSTR_Exp_ | pcGEM integrated with experimental proteomics for suboptimal condition, used for constructing the baseline and calculating correlations and RMdSE |
| _lambda_      | weight for the sum of the distance between enzyme allocations and the distance between flux distributions |
