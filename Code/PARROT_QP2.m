%
% Minimizes ||Eref-Es||2 + λ||Vref-Vs||2
%
% Code partially adapted from MOMA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
changeCobraSolver('gurobi', 'QP');

%% Loading the enzyme-constrained model and other data
load('modelREF.mat')                    % pcGEM integrated with experimental proteomics for reference condition 
load('modelSTR.mat')                    % pcGEM without integrated proteomics data
load("modelSTR_Exp.mat")                % pcGEM integrated with experimental proteomics for suboptimal condition
load('baseline_filename.mat')           % Baseline for calculating comparisons

%% Constraints for generating Es
modelSTR = changeRxnBounds(modelSTR, 'glucose', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'CO2', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'O2', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'pyruvate', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'succinate', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'glycerol', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'acetate', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'ethanol', 1*(1+0.5), 'u');

modelSTR = changeRxnBounds(modelSTR, 'glutamate', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'phenylalanine', 1*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'isoleucine', 1*(1+0.5), 'u');

% fix specific growth rate at the dilution rate, allowing 5% flexibility
modelSTR = changeRxnBounds(modelSTR, 'biomass reaction', 0.1*(1+0.05),'u');
modelSTR = changeRxnBounds(modelSTR, 'biomass reaction', 0.1*(1-0.05),'l');

ptotREF = 0.5;
ptotSTR = 0.5;
lambda = 1;

%% Generating VS1 and ES1 (REF)
[~,nRxns] = size(modelREF.S);
enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
enzymeIds(end) = [];

infNum = sum(isinf(modelREF.ub(enzymeIds)));
infBounds = find(isinf(modelREF.ub));
infBoundsProt = intersect(infBounds, enzymeIds);

ptotShare = modelREF.ub(nRxns)/infNum;
modelREF.ub(infBoundsProt) = ptotShare;

NormMin = buildLPproblemFromModel(modelREF);

proteins = modelREF.enzymes(:);

pIdx = [];
rxnIdx = {};
kcat = [];

for i=1:numel(proteins)
    pIdx(i,1) = find(strcmpi(modelREF.metNames,join(['prot_' char(proteins(i))],"")));
    rxnIdx{i,1} = find(modelREF.S(pIdx(i),:) < 0);
    kcat(i,1) = min(-1./modelREF.S(pIdx(i),rxnIdx{i}));
end

NormMin.c = zeros(size(NormMin.A,2),1);
NormMin.c(enzymeIds) = 1;
NormMin.c(enzymeIds) = kcat.*NormMin.c(enzymeIds);

NormMin.osense = 1;

fprintf('\n');
validLPProblem = verifyCobraProblem(NormMin);

NormMinSol = solveCobraLP(NormMin);

%% Construct the QP problem to find Es
[~,nRxnsSTR] = size(modelSTR.S);
enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
[enzymeNum,~] = size(enzymeIds);
rxnIds = find(~cellfun('isempty',modelREF.rxns(1:min(enzymeIds)-1)));

ptotSTR = ptotREF;

QPproblem = buildLPproblemFromModel(modelSTR);

QPproblem.c(rxnIds) = NormMinSol.full(rxnIds);
QPproblem.c(rxnIds) = -2*lambda*QPproblem.c(rxnIds);
QPproblem.c(enzymeIds) = NormMinSol.full(enzymeIds)/ptotREF;
QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);

QPproblem.F = sparse(size(QPproblem.A,2));
QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 2*lambda*speye(nRxnsSTR);
QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);

QPproblem.osense = 1;

%% Verify if the problem is a valid QP problem
fprintf('\n');
validQPProblem = verifyCobraProblem(QPproblem);

%% Solve the problem (COBRA)
QPsolution = solveCobraQP(QPproblem);

%% Get the solution(s)
MinimizedFlux.x = QPsolution.full;
  
MinimizedFlux.f = sum(modelSTR.c.*MinimizedFlux.x);

MinimizedFlux.stat = QPsolution.stat;
solStatus = QPsolution.stat;
MinimizedFlux.solver = QPsolution.solver;
MinimizedFlux.time = QPsolution.time;

MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds)*ptotSTR;

fprintf('\n');
fprintf('QP SOLUTION\n');
printFluxes(modelSTR, MinimizedFlux.x, true);

fprintf('\n');
fprintf('QP SOLUTION FOR ES2 \n');
protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
protNames.names = modelSTR.metNames(protNames.ids);
printFluxes(modelSTR, MinimizedFlux.x, false, '', '', '', protNames.names);

%% Calculate correlation with experimental values
enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];

pred = {};
pred(:,1) = modelSTR.rxns(enzymeIds);
pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
pred(:,2) = num2cell(MinimizedFlux.x(enzymeIds));
pred = cell2table(pred);
pred.Properties.VariableNames = {'Protein' 'Predicted'};
pred.Protein = char(pred.Protein);

merged = innerjoin(pred, baseline);
merged(ismember(merged.Predicted, 0),:)=[];
merged(ismember(merged.Abundance, 0),:)=[];

merged.Predicted = log10(abs(merged.Predicted));
merged.Abundance = log10(abs(merged.Abundance));

cor = {};
cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};

[cor{2,1}, cor{3,1}] = corr(merged.Abundance, merged.Predicted, 'Type', char(cor(1,1)));
[cor{2,2}, cor{3,2}] = corr(merged.Abundance, merged.Predicted, 'Type', char(cor(1,2)));
[cor{2,3}, cor{3,3}] = corr(merged.Abundance, merged.Predicted, 'Type', char(cor(1,3)));

formatSpec = "%s's correlation: %f \t p-value: %s \n";
fprintf('\n');
fprintf('Correlation between predicted and experimental values:')
fprintf('\n');
fprintf(formatSpec, cor{:});

%% Calculate the RMSE
rmse = cell2table(cell(0,0));
rmse.Protein = modelSTR_Exp.enzymes;
rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
rmse.Predicted = log10(abs(MinimizedFlux.x(enzymeIds)));

rmse.Experimental(isinf(rmse.Experimental)) = NaN;
rmse.Predicted(isinf(rmse.Predicted)) = NaN;
rmse = rmmissing(rmse);

rmse.RMSE = rmse.Experimental - rmse.Predicted;
rmse.RMSE = (rmse.RMSE).^2;
result = median(rmse.RMSE);
result = sqrt(result);

fprintf('\n');
formatSpecRMSE = "Log10 root median squared error: %f";
fprintf(formatSpecRMSE, result);
fprintf('\n');

%% How many proteins inside Es?
fprintf('\n');
[numEs,~] = size(merged);
formatSpecNUM = "Number of predicted Es: %f";
fprintf(formatSpecNUM, numEs);
fprintf('\n');

%% Plot scatterplot of predicted values vs. baseline
figure
hold on
scatter(merged.Abundance,merged.Predicted, 'filled')
title(sp_title)
xlabel('Experimental values')
ylabel('Predicted values')
p = polyfit(merged.Abundance, merged.Predicted, 1);
px = [min(merged.Abundance) max(merged.Abundance)];
py = polyval(p, px);
R = corrcoef(merged.Abundance, merged.Predicted);
str = join(['r = ', char(num2str(R(2)))], "");
text(-5.5, -11.5, str);
plot(px, py, 'LineWidth', 2);
hold off