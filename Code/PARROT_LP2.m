%
% Minimizes ||Eref-Es|| + λ||Vref-Vs||
%
% Code partially adapted from linearMOMA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
changeCobraSolver('gurobi', 'LP');

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
[enzymeNum,~] = size(enzymeIds);

infNum = sum(isinf(modelREF.ub(enzymeIds)));
infBounds = find(isinf(modelREF.ub));
infBoundsProt = intersect(infBounds, enzymeIds);

ptotShare = (modelREF.ub(nRxns)/infNum);
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

% printFluxes(modelREF, NormMinSol.full, true);

%% Construct the LP problem to find ES2
[~,nRxnsREF] = size(modelREF.S);
[~,nRxnsSTR] = size(modelSTR.S);
enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
[enzymeNum,~] = size(enzymeIds);

LPREF = buildLPproblemFromModel(modelREF);
LPSTR = buildLPproblemFromModel(modelSTR);

[nREFCtrs,nREFVars] = size(LPREF.A);
[nSTRCtrs,nSTRVars] = size(LPSTR.A);

commonRxns = ismember(modelREF.rxns,modelSTR.rxns);
nCommon = sum(commonRxns);

deltaMat = createDeltaMatchMatrix(modelREF.rxns,modelSTR.rxns);
deltaMatREF = deltaMat(1:2*nCommon,1:nRxnsREF);
deltaMatSTR = deltaMat(1:2*nCommon,nRxnsREF+(1:nRxnsSTR));
deltaMatCom = deltaMat(1:2*nCommon,(nRxnsREF+nRxnsSTR)+(1:2*nCommon));

LPproblem.A = [sparse(nSTRCtrs,nREFVars),LPSTR.A,sparse(nSTRCtrs,2*nCommon);...
    deltaMatREF, sparse(2*nCommon,nREFVars - nRxnsREF), deltaMatSTR, sparse(2*nCommon,nSTRVars - nRxnsSTR), deltaMatCom;...
    LPREF.c',sparse(1,nSTRVars+2*nCommon)];

% Construct the RHS vector
LPproblem.b = [LPSTR.b;zeros(2*nCommon,1);1];

% Linear objective = 0
enzymeIds_cindex = enzymeIds + nREFVars;
rxns_cindex = nSTRVars - enzymeNum;

LPproblem.c = [zeros(nREFVars+nSTRVars,1); zeros(2*nCommon,1)];
LPproblem.c(1:rxns_cindex) = lambda*NormMinSol.full(1:rxns_cindex);
LPproblem.c(enzymeIds_cindex) = (NormMinSol.full(enzymeIds)/ptotREF);

% Construct the ub/lb
% delta [-10000 10000]
LPproblem.lb = [zeros(nREFVars,1);LPSTR.lb;zeros(2*nCommon,1)];
LPproblem.ub = [zeros(nREFVars,1);LPSTR.ub;10000*ones(2*nCommon,1)];

% minimize
LPproblem.osense = 1;

% Construct the constraint direction vector (G for delta's, E for
% everything else)
csense = 'L';

LPproblem.csense = [LPSTR.csense;repmat('G',2*nCommon,1);csense];

%% Verify if the problem is a valid QP problem
fprintf('\n');
statusOK = verifyCobraProblem(LPproblem);

%% Solve the problem
MinimizedFlux = solveCobraLP(LPproblem);
MinimizedFlux.x = MinimizedFlux.full(8144:16287);
MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * ptotSTR;

%% Get the solution(s)
fprintf('\n');
fprintf('LP SOLUTION\n');
printFluxes(modelSTR, MinimizedFlux.x, true);

fprintf('\n');
fprintf('LP SOLUTION FOR ES2 \n');
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