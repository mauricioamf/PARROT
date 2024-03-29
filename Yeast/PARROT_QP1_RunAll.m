%
% Minimizes ||Eref-Es||2 + λ||Vref-Vs||2
%
% Code partially adapted from MOMA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
tic
changeCobraSolver('gurobi', 'QP');

%% Loading the enzyme-constrained model and other data
% Constraints for suboptimal conditions
load('./constraints_Scerevisiae.mat')

% We need to calculate the total enzyme usage Etot, which is the total
% protein content Ptot, multiplied by the fraction f of enzymes that are 
% accounted for in the model, and a parameter sigma, which is the average 
% in vivo saturation of all enzymes. This is done for both reference and
% alternative models. For both models, Etot is used to normalize the enzyme
% usage distributions. For the alternative growth model, Etot is used to
% constrain the protein pool exchange pseudo-reaction as well.
%
% The values of sigma and f can be specified, but in their absence a
% default of 0.5 is usually assumed. We noticed that 0.5 resulted in 
% unfeasible solutions. We adopted 0.4 as the default because this was 
% the value that allowed growth while ensuring that the protein pool 
% exchange pseudo-reaction was adequately constrained.
%
% For all the other parameters used in this code (uptake rates, growth 
% rates, Ptots, etc), refer to the file "constraints_Scerevisiae.xlsx".
%
% In case a particular growth condition results in unfeasible solutions, we
% suggest running PARROT for that individual condition, using the function
% PARROT.m and constraining the model as shown in the examples. The
% constraints can be found in the file "constraints_Scerevisiae.xlsx".

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
growth_conditions = constraints.Properties.VariableNames(3:9);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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

ptotREF = 0.265;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

results_YU2021_N30 = cell2table(cell(0,0));
% corrVals_YU2021_N30 = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(10:15);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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
growth_conditions = constraints.Properties.VariableNames(17:19);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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
growth_conditions = constraints.Properties.VariableNames(20:21);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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
growth_conditions = constraints.Properties.VariableNames(22:23);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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
growth_conditions = constraints.Properties.VariableNames(24:27);

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
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
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

    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelSTR);
    
    QPproblem.c = zeros(nRxnsSTR,1);
    QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
    QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
    QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
    
    QPproblem.F = sparse(size(QPproblem.A,2));
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = speye(nRxnsSTR);
    QPproblem.F(1:nRxnsSTR,1:nRxnsSTR) = 0;
    QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
    
    QPproblem.osense = 1;

    % Solve the problem
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
    MinimizedFlux.x = QPsolution.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds) * etotSTR;
    
    % fprintf('\n');
    % fprintf('QP SOLUTION\n');
    % printFluxes(modelSTR, MinimizedFlux.x, true);
    
    % fprintf('\n');
    % fprintf('QP SOLUTION FOR ES2 \n');
    % protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    % protNames.names = modelSTR.metNames(protNames.ids);
    % filename = current_condition + ".txt";
    % printFluxes(modelSTR, MinimizedFlux.x, false, '', filename, '', protNames.names);
    
    % Calculate correlation with experimental values
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
    % fprintf("Correlation between predicted and experimental values:")
    % fprintf('\n');
    % fprintf(formatSpec, cor{:});
       
    % Calculate the RMSE
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
    
    % fprintf('\n');
    % formatSpecRMSE = "Log10 root median squared error: %f";
    % fprintf(formatSpecRMSE, result);
    % fprintf('\n');
    
    % How many proteins inside Es?
    fprintf('\n');
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
    % text(-5.5, -11.5, str);
    % plot(px, py, 'LineWidth', 2);
    % title(current_condition);
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

% corrVals_final = [corrVals_LAHTVEE;
%     corrVals_YU2021_N30;
%     corrVals_YU2021_Gln;
%     corrVals_YU2021_Phe;
%     corrVals_YU2021_Ile;
%     corrVals_YU2020];

%%
toc