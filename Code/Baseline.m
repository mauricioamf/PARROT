%
% Constructs baseline by minimizing the first norm of the enzyme usage
% distribution
%
%%
clear;clc
changeCobraSolver('gurobi', 'LP');

%%
load("modelSTR_Exp.mat") % pcGEM integrated with experimental proteomics for suboptimal condition

baseline_filename = 'baseline_filename.mat'; % filename to export the baseline file

%%
% fix specific growth rate at the dilution rate, allowing 5% flexibility
model = changeRxnBounds(model, 'biomass reaction', 0.1*(1+0.05),'u');
model = changeRxnBounds(model, 'biomass reaction', 0.1*(1-0.05),'l');
 
% substrate uptake, allowing 5% flexibility
model = changeRxnBounds(model, 'substrate_REV', 0.1*(1+0.5),'u');   % constrain the REV reaction
model = changeRxnBounds(model, 'substrate_REV', 0.1*(1-0.5),'l');
model = changeRxnBounds(model, 'substrate', 0,'b'); % block the original reaction

% additional constraints, allowing 5% flexibility
model = changeRxnBounds(model, 'glucose', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'CO2', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'O2', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'pyruvate', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'succinate', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'glycerol', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'acetate', 1*(1+0.5), 'u');
model = changeRxnBounds(model, 'ethanol', 1*(1+0.5), 'u');

% minimize enzyme usage
model = changeObjective(model, 'biomass reaction', 0);
enzyme_usage = find(~cellfun('isempty',strfind(model.rxnNames,'prot_')));
model = changeObjective(model, model.rxns(enzyme_usage), 1);

%%
solBL = optimizeCbModel(model, 'min');

%%
fprintf('\n');
printFluxes(model, solBL.x, true);

%%
enzymeIds = find(~cellfun('isempty',strfind(model.rxnNames,'prot_'))); 

baseline = {};
baseline(:,1) = model.rxns(enzymeIds);
baseline(:,1) = replace(baseline(:,1), 'draw_prot_', '');
baseline(:,1) = replace(baseline(:,1), 'prot_', '');
baseline(:,1) = replace(baseline(:,1), '_exchange', '');
baseline(:,2) = num2cell(solBL.x(enzymeIds));
baseline = cell2table(baseline);
baseline.Properties.VariableNames = {'Protein' 'Abundance'};
baseline.Protein = char(baseline.Protein);
baseline(end,:) = [];

save(baseline_filename, 'baseline');