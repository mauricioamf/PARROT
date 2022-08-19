%
% Constructs baseline by minimizing the second norm of the enzyme usage
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

%%
[~,nRxns] = size(model.S);
enzymeIds = find(~cellfun('isempty',strfind(model.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];
[enzymeNum,~] = size(enzymeIds);

QPproblem = buildLPproblemFromModel(model);

QPproblem.c = zeros(size(QPproblem.A,2),1);

QPproblem.F = sparse(size(QPproblem.A,2));
QPproblem.F(1:nRxns,1:nRxns) = 2*speye(nRxns);
QPproblem.F(1:nRxns,1:nRxns) = 0;
QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);

QPproblem.osense = 1;

%% Verify if the problem is a valid QP problem
fprintf('\n');
validQPProblem = verifyCobraProblem(QPproblem);

%% Solve the problem (COBRA)
QPsolution = solveCobraQP(QPproblem);

%% Get the solution(s)
MinimizedFlux.x = QPsolution.full;
  
MinimizedFlux.f = sum(model.c.*MinimizedFlux.x);

solStatus = QPsolution.stat;
MinimizedFlux.stat = QPsolution.stat;
MinimizedFlux.solver = QPsolution.solver;
MinimizedFlux.time = QPsolution.time;

fprintf('\n');
fprintf('QP SOLUTION\n');
printFluxes(model, MinimizedFlux.x, true);

%%
enzymeIds = find(~cellfun('isempty',strfind(model.rxnNames,'prot_'))); 

baseline = {};
baseline(:,1) = model.rxns(enzymeIds);
baseline(:,1) = replace(baseline(:,1), 'draw_prot_', '');
baseline(:,1) = replace(baseline(:,1), 'prot_', '');
baseline(:,1) = replace(baseline(:,1), '_exchange', '');
baseline(:,2) = num2cell(MinimizedFlux.x(enzymeIds));
baseline = cell2table(baseline);
baseline.Properties.VariableNames = {'Protein' 'Abundance'};
baseline.Protein = char(baseline.Protein);
baseline(end,:) = [];

save(baseline_filename, 'baseline');