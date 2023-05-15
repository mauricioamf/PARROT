%
% Performs pFBA
%
% Code partially adapted from pFBA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
tic
changeCobraSolver('ibm_cplex', 'LP');

%% Loading the enzyme-constrained model and other data
% Constraints for suboptimal conditions
load('./constraints_Scerevisiae_noref.mat')

%% Loop over growth conditions using Lahtvee2017_REF as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Lahtvee2017_REF.mat')                  
modelREF = model;
clear model

ptotREF = 0.4228;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_LAHTVEE = cell2table(cell(0,0));
% corrVals_LAHTVEE = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:8);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.05);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_LAHTVEE.Conditions(k) = growth_conditions(k);
    results_LAHTVEE.Pearson(k) = cor{2,1};
    results_LAHTVEE.pvalue_P(k) = cor{3,1};
    results_LAHTVEE.Spearman(k) = cor{2,2};
    results_LAHTVEE.pvalue_S(k) = cor{3,2};
    results_LAHTVEE.RMdSE(k) = result;
    results_LAHTVEE.NumEs(k) = numEs;

end

%% Loop over growth conditions using Yu2021_std_010 as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Yu2021_std_010.mat')                  
modelREF = model;
clear model

ptotREF = 0.3665;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2021_N30 = cell2table(cell(0,0));
% corrVals_YU2021_N30 = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(9:13);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.5);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.5), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_YU2021_N30.Conditions(k) = growth_conditions(k);
    results_YU2021_N30.Pearson(k) = cor{2,1};
    results_YU2021_N30.pvalue_P(k) = cor{3,1};
    results_YU2021_N30.Spearman(k) = cor{2,2};
    results_YU2021_N30.pvalue_S(k) = cor{3,2};
    results_YU2021_N30.RMdSE(k) = result;
    results_YU2021_N30.NumEs(k) = numEs;

end

%% Loop over growth conditions using Yu2021_Gln_glc1 as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Yu2021_Gln_glc1.mat')                  
modelREF = model;
clear model

ptotREF = 0.3665;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2021_Gln = cell2table(cell(0,0));
% corrVals_YU2021_Gln = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(15:16);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.5);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.5), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_YU2021_Gln.Conditions(k) = growth_conditions(k);
    results_YU2021_Gln.Pearson(k) = cor{2,1};
    results_YU2021_Gln.pvalue_P(k) = cor{3,1};
    results_YU2021_Gln.Spearman(k) = cor{2,2};
    results_YU2021_Gln.pvalue_S(k) = cor{3,2};
    results_YU2021_Gln.RMdSE(k) = result;
    results_YU2021_Gln.NumEs(k) = numEs;

end

%% Loop over growth conditions using Yu2021_Phe_std as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Yu2021_Phe_std.mat')                  
modelREF = model;
clear model

ptotREF = 0.5586;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2021_Phe = cell2table(cell(0,0));
% corrVals_YU2021_Phe = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(17:17);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.5);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.5), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_YU2021_Phe.Conditions(k) = growth_conditions(k);
    results_YU2021_Phe.Pearson(k) = cor{2,1};
    results_YU2021_Phe.pvalue_P(k) = cor{3,1};
    results_YU2021_Phe.Spearman(k) = cor{2,2};
    results_YU2021_Phe.pvalue_S(k) = cor{3,2};
    results_YU2021_Phe.RMdSE(k) = result;
    results_YU2021_Phe.NumEs(k) = numEs;

end

%% Loop over growth conditions using Yu2021_Ile_std as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Yu2021_Ile_std.mat')                  
modelREF = model;
clear model

ptotREF = 0.49431;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2021_Ile = cell2table(cell(0,0));
% corrVals_YU2021_Ile = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(18:18);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.5);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.5), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_YU2021_Ile.Conditions(k) = growth_conditions(k);
    results_YU2021_Ile.Pearson(k) = cor{2,1};
    results_YU2021_Ile.pvalue_P(k) = cor{3,1};
    results_YU2021_Ile.Spearman(k) = cor{2,2};
    results_YU2021_Ile.pvalue_S(k) = cor{3,2};
    results_YU2021_Ile.RMdSE(k) = result;
    results_YU2021_Ile.NumEs(k) = numEs;

end

%% Loop over growth conditions using Yu2020_Clim as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/Yu2020_Clim.mat')                  
modelREF = model;
clear model

ptotREF = 0.4658;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2020 = cell2table(cell(0,0));
% corrVals_YU2020 = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(19:21);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/ecYeastGEM_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = model;
    clear model

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir)  
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.05);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(13);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Construct the LP problem to find Es
    [~,nRxnsSTR] = size(modelSTR.S);
    gprAssociated = find(sum(modelSTR.rxnGeneMat,2)>0);
    growth = find(contains(modelSTR.rxnNames, 'biomass pseudoreaction'));
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    
    modelSTR = convertToIrreversible(modelSTR);
    
    LPproblem = buildLPproblemFromModel(modelSTR);
    
    LPproblem.lb(LPproblem.lb==-Inf) = -1000;
    LPproblem.ub(LPproblem.ub==Inf) = 1000;
    
    enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    
    LPproblem.c = zeros(nRxnsSTR,1);
    LPproblem.c(gprAssociated) = 1;
    %LPproblem.c(min(enzymeIds):max(enzymeIds)-1) = modelREF.ub(min(enzymeIds):max(enzymeIds)-1);
    
    LPproblem.osense = 1;
    
    % Solve the problem
    LPsolution = solveCobraLP(LPproblem);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full;
    LPsolution.x(enzymeIds) = LPsolution.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(modelSTR, LPsolution.x, true);
    
    % fprintf('\n');
    % fprintf('LP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % printFluxes(modelSTR, LPsolution.x, false, '', '', '', protNames.names);
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(LPsolution.x(enzymeIds));
    pred = cell2table(pred);
    pred.Properties.VariableNames = {'Protein' 'Predicted'};
    pred.Protein = char(pred.Protein);
    
    merged = innerjoin(pred, baseline);
    merged(ismember(merged.Predicted, 0),:)=[];
    merged(ismember(merged.Abundance, 0),:)=[];
    
    mergedLog = merged;
    mergedLog.Predicted = log10(abs(mergedLog.Predicted));
    mergedLog.Abundance = log10(abs(mergedLog.Abundance));
    
    cor = {};
    cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};
    
    [cor{2,1}, cor{3,1}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,2)));
    [cor{2,3}, cor{3,3}] = corr(mergedLog.Abundance, mergedLog.Predicted, 'Type', char(cor(1,3)));
    
    % formatSpec = "%s's correlation: %f \t p-value: %s \n";
    % fprintf('\n');
    % fprintf('Correlation between predicted and experimental values:')
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
    
    % Calculate the RMSE
    rmse = cell2table(cell(0,0));
    rmse.Protein = modelSTR_Exp.enzymes;
    rmse.Experimental = log10(abs(modelSTR_Exp.ub(enzymeIds)));
    rmse.Predicted = log10(abs(LPsolution.x(enzymeIds)));
    
    rmse.Experimental(isinf(rmse.Experimental)) = NaN;
    rmse.Predicted(isinf(rmse.Predicted)) = NaN;
    rmse = rmmissing(rmse);
    
    rmse.RMSE = rmse.Experimental - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    result = median(rmse.RMSE);
    result = sqrt(result);
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    % fprintf('\n');
    [numEs,~] = size(merged);
    % formatSpecNUM = "Number of predicted Es: %f";
    % fprintf(formatSpecNUM, numEs);
    % fprintf('\n');
    
    % Plot scatterplot of predicted values vs. baseline
    % figure
    % hold on
    % scatter(mergedLog.Abundance,mergedLog.Predicted, 'filled')
    % xlabel('Experimental values')
    % ylabel('Predicted values')
    % p = polyfit(mergedLog.Abundance, mergedLog.Predicted, 1);
    % px = [min(mergedLog.Abundance) max(mergedLog.Abundance)];
    % py = polyval(p, px);
    % R = corrcoef(mergedLog.Abundance, mergedLog.Predicted);
    % str = join(['r = ', char(num2str(R(2)))], "");
    % text(-5, -8, str);
    % plot(px, py, 'LineWidth', 2);
    % hold off

    results_YU2020.Conditions(k) = growth_conditions(k);
    results_YU2020.Pearson(k) = cor{2,1};
    results_YU2020.pvalue_P(k) = cor{3,1};
    results_YU2020.Spearman(k) = cor{2,2};
    results_YU2020.pvalue_S(k) = cor{3,2};
    results_YU2020.RMdSE(k) = result;
    results_YU2020.NumEs(k) = numEs;

end

%% Merge all results tables
results_table = [results_LAHTVEE; 
    results_YU2021_N30; 
    results_YU2021_Gln; 
    results_YU2021_Phe;
    results_YU2021_Ile;
    results_YU2020];

%%
toc