%% EXAMPLE 1 - How does S. cerevisiae reallocates enzyme usage under osmotic stress?
% This file contains an example that demonstrates the functionality of
% PARROT. The example shows how to run PARROT using default
% parameters to predict the enzyme usage distribution of S. cerevisae under
% high osmolarity stress. This example uses the ecYeast8 model (DOI: 10.1038/s41467-019-11581-3)
% and proteomics data from Lahtvee et al. (2017, DOI: 10.1016/j.cels.2017.03.003).
% In this example, we use REF to refer to the reference growth condition,
% and ALT to refer to the alternative growth condition, and variables have
% been named accordingly.
%
% The function PARROT.m was designed to work similarly to the MOMA.m
% function included in the COBRA Toolbox.
%
% Author: Mauricio Ferreira, 2023

%% Preparation before running the examples

% Cleaning the workspace and the command window
clear;clc;

% Initialize COBRA and set the solver to Gurobi
% initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
changeCobraSolver('gurobi', 'QP');

%% STEP 1: Load the enzyme-constrained models

% First, load the pcGEM integrated with experimental proteomics for 
% the reference condition. The enzyme usage distribution of this model is
% used as the starting point for predicting the enzyme usage distribution
% of the alternative condition.
load('./ExampleData/Lahtvee2017_REF.mat')                  
modelREF = model;
clear model

% Then, load the batch pcGEM, which is the version without any integrated 
% proteomics data. This model is used to predict the enzyme usage under
% the alternative growth condition. 
load('./ExampleData/ecYeastGEM_batch.mat')          
modelALT = ecModel_batch;
clear ecModel_batch

%% STEP 2: Apply constraints 

% Here we apply the constraints that represents the alternative growth
% condition. The predictions work better when experimental data is used in
% this step. If the model ends up overconstrained, then relaxing the
% constraints is required
modelALT = changeRxnBounds(modelALT, 'r_1714_REV', 1.36797731358331, 'u');   % Glucose
modelALT = changeRxnBounds(modelALT, 'r_1672', 4.11867992630463, 'u');       % CO2
modelALT = changeRxnBounds(modelALT, 'r_1992_REV', 4.08889827981797, 'u');   % O2
modelALT = changeRxnBounds(modelALT, 'r_1808', 0.0283248589619358, 'u');     % Glycerol

% Fix specific growth rate at the dilution rate
modelALT = changeRxnBounds(modelALT, 'r_2111', 0.1, 'u');
modelALT = changeRxnBounds(modelALT, 'r_2111', 0.1, 'l');

% Next, need to calculate the total enzyme usage Etot, which is the total
% protein content Ptot, multiplied by the fraction f of enzymes that are 
% accounted for in the model, and a parameter sigma, which is the average 
% in vivo saturation of all enzymes. This is done for both reference and
% alternative models. For both models, Etot is used to normalize the enzyme
% usage distributions. For the alternative growth model, Etot is used to
% constrain the protein pool exchange pseudo-reaction as well.
ptotREF = 0.4228; % Lahtvee et al (2017, DOI: 10.1016/j.cels.2017.03.003) 
ptotALT = 0.64;   % Lahtvee et al (2017, DOI: 10.1016/j.cels.2017.03.003)

% The values of sigma and f can be specified, but in their absence a
% default of 0.5 can be assumed.
sigmaREF = 0.5;
fREF = 0.5;
sigmaALT = 0.5;
fALT = 0.5;

% Finally, we can get Etot for both conditions
etotREF = ptotREF * sigmaREF * fREF;
etotALT = ptotALT * sigmaALT * fALT;

% And then constrain the protein pool exchange pseudo-reaction for the
% alternative condition
modelALT = changeRxnBounds(modelALT, 'prot_pool_exchange', etotALT, 'b');

%% STEP 3: Run PARROT

% Now, all we need to do is run the function. The input parameters for
% PARROT are:
%
% - modelREF: the reference growth condition model
%
% - modelALT: the alternative growth condition model
%
% - MinStr: the minimization strategy. In this example, we use the default
% 'Manhattan'
%
% - lambda: the weighting factor for metabolic fluxes. In this example, we
% don't use fluxes, only the enzyme usage distributions
%
% - etotREF: the total enzyme usage for the reference growth condition
% model. This parameter is used to normalize the enzyme usage distribution
%
% - etotALT: the total enzyme usage for the alternative growth condition
% model. This parameter normalizes the enzyme usage distribution and also
% constrains the protein pool exchange pseudo-reaction
%
% In this example, the parameters are set to represent the LP1 variant
% described in the manuscript, which minimizes the Manhattan distance of
% the enzyme usage distribution of the reference and alternative growth
% conditions, and don't consider metabolic fluxes in the problem

[solutionALT, solStatus] = PARROT(modelREF, modelALT, 'Manhattan', 0, etotREF, etotALT);

%% STEP 4: Get the solution(s)

% Let's check the flux distribution for the alternative growth condition:
fprintf('\n');
fprintf('PARROT SOLUTION FOR FLUXES \n');
printFluxes(modelALT, solutionALT, true);

% Next, let's check the enzyme usage distribution for the alternative
% growth condition:
fprintf('\n');
fprintf('PARROT SOLUTION FOR ES \n');
protNames.ids = find(~cellfun('isempty',strfind(modelALT.metNames,'prot_pool')));
protNames.names = modelALT.metNames(protNames.ids);
printFluxes(modelALT, solutionALT, false, '', '', '', protNames.names);