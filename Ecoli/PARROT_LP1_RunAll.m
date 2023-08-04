%
% Minimizes ||Eref-Es||
%
% Code partially adapted from MOMA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
tic
changeCobraSolver('gurobi', 'LP');

%% Loading the enzyme-constrained model and other data
% Constraints for suboptimal conditions
load('./constraints_Ecoli.mat')

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
% rates, Ptots, etc), refer to the file "constraints_Ecoli.xlsx".
%
% In case a particular growth condition results in unfeasible solutions, we
% suggest running PARROT for that individual condition, using the function
% PARROT.m and constraining the model as shown in the examples. The
% constraints can be found in the file "constraints_Ecoli.xlsx".

ptotREF = 0.61;

f_REF = 0.4;
sigma_REF = 0.4;
f_STR = 0.4;
sigma_STR = 0.4;

%% Loop over growth conditions using GLC_BATCH_mu_0_58_S as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/GLC_BATCH_mu_0_58_S.mat')                  
modelREF = ecModelP;
clear ecModelP

results_BATCH = cell2table(cell(0,0));
% corrVals_CHEM_S = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:9);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/eciML1515_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
    load(baselinedir)  
    baseline = Ecoli;
    clear Ecoli
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(9);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % Construct the LP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
    [enzymeNum,~] = size(enzymeIds);
    
    % infNum = sum(isinf(modelREF.ub(enzymeIds)));
    % infBounds = find(isinf(modelREF.ub));
    % infBoundsProt = intersect(infBounds, enzymeIds);
    % 
    % etotShare = etotSTR/infNum;
    % modelREF.ub(infBoundsProt) = etotShare;
    
    modelDelta = struct();
    
    % Set the lower and upper bounds for the fluxes
    modelDelta.lb = [modelREF.lb; modelSTR.lb];
    modelDelta.lb(modelDelta.lb==-Inf) = -1000;
    
    modelDelta.ub = [modelREF.ub; modelSTR.ub];
    modelDelta.ub(modelDelta.ub==Inf) = 1000;
    
    % Concatenate the stoichiometric matrices and S-matrices of the two models
    modelDelta.S = [modelREF.S zeros(size(modelREF.S,1), size(modelSTR.S,2)); zeros(size(modelSTR.S,1), size(modelREF.S,2)) modelSTR.S];
    
    % Set the RHS vector
    modelDelta.b = zeros(size(modelDelta.S,1),1);
    
    % Objective function
    [~,nVars] = size(modelDelta.S);
    % enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    enzymeIds_cindex = enzymeIds + nRxnsREF;
    
    modelDelta.c = zeros(nVars,1);
    modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
    modelDelta.c(enzymeIds) = modelDelta.c(enzymeIds)/etotREF;
    modelDelta.c(enzymeIds_cindex) = -ones(length(modelSTR.rxns(enzymeIds)),1);
    
    % Solve the problem
    LPsolution = solveCobraLP(modelDelta);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full(length(modelREF.rxns)+1:end);
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

    results_BATCH.Conditions(k) = growth_conditions(k);
    results_BATCH.Pearson(k) = cor{2,1};
    results_BATCH.pvalue_P(k) = cor{3,1};
    results_BATCH.Spearman(k) = cor{2,2};
    results_BATCH.pvalue_S(k) = cor{3,2};
    results_BATCH.RMdSE(k) = result;
    results_BATCH.NumEs(k) = numEs;

end

%% Loop over growth conditions using GLC_CHEM_mu_0_12_S as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/GLC_CHEM_mu_0_12_S.mat')                  
modelREF = ecModelP;
clear ecModelP

results_CHEM_S = cell2table(cell(0,0));
% corrVals_CHEM_S = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(10:13);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/eciML1515_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
    load(baselinedir)  
    baseline = Ecoli;
    clear Ecoli
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(9);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % Construct the LP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
    [enzymeNum,~] = size(enzymeIds);
    
    % infNum = sum(isinf(modelREF.ub(enzymeIds)));
    % infBounds = find(isinf(modelREF.ub));
    % infBoundsProt = intersect(infBounds, enzymeIds);
    % 
    % etotShare = etotSTR/infNum;
    % modelREF.ub(infBoundsProt) = etotShare;
    
    modelDelta = struct();
    
    % Set the lower and upper bounds for the fluxes
    modelDelta.lb = [modelREF.lb; modelSTR.lb];
    modelDelta.lb(modelDelta.lb==-Inf) = -1000;
    
    modelDelta.ub = [modelREF.ub; modelSTR.ub];
    modelDelta.ub(modelDelta.ub==Inf) = 1000;
    
    % Concatenate the stoichiometric matrices and S-matrices of the two models
    modelDelta.S = [modelREF.S zeros(size(modelREF.S,1), size(modelSTR.S,2)); zeros(size(modelSTR.S,1), size(modelREF.S,2)) modelSTR.S];
    
    % Set the RHS vector
    modelDelta.b = zeros(size(modelDelta.S,1),1);
    
    % Objective function
    [~,nVars] = size(modelDelta.S);
    % enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    enzymeIds_cindex = enzymeIds + nRxnsREF;
    
    modelDelta.c = zeros(nVars,1);
    modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
    modelDelta.c(enzymeIds) = modelDelta.c(enzymeIds)/etotREF;
    modelDelta.c(enzymeIds_cindex) = -ones(length(modelSTR.rxns(enzymeIds)),1);
    
    % Solve the problem
    LPsolution = solveCobraLP(modelDelta);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full(length(modelREF.rxns)+1:end);
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

    results_CHEM_S.Conditions(k) = growth_conditions(k);
    results_CHEM_S.Pearson(k) = cor{2,1};
    results_CHEM_S.pvalue_P(k) = cor{3,1};
    results_CHEM_S.Spearman(k) = cor{2,2};
    results_CHEM_S.pvalue_S(k) = cor{3,2};
    results_CHEM_S.RMdSE(k) = result;
    results_CHEM_S.NumEs(k) = numEs;

end

%% Loop over growth conditions using GLC_CHEM_mu_0_11_V as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/GLC_CHEM_mu_0_11_V.mat')                  
modelREF = ecModelP;
clear ecModelP

results_CHEM_V = cell2table(cell(0,0));
% corrVals_CHEM_V = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(14:18);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/eciML1515_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
    load(baselinedir)  
    baseline = Ecoli;
    clear Ecoli
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(9);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % Construct the LP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
    [enzymeNum,~] = size(enzymeIds);
    
    % infNum = sum(isinf(modelREF.ub(enzymeIds)));
    % infBounds = find(isinf(modelREF.ub));
    % infBoundsProt = intersect(infBounds, enzymeIds);
    % 
    % etotShare = etotSTR/infNum;
    % modelREF.ub(infBoundsProt) = etotShare;
    
    modelDelta = struct();
    
    % Set the lower and upper bounds for the fluxes
    modelDelta.lb = [modelREF.lb; modelSTR.lb];
    modelDelta.lb(modelDelta.lb==-Inf) = -1000;
    
    modelDelta.ub = [modelREF.ub; modelSTR.ub];
    modelDelta.ub(modelDelta.ub==Inf) = 1000;
    
    % Concatenate the stoichiometric matrices and S-matrices of the two models
    modelDelta.S = [modelREF.S zeros(size(modelREF.S,1), size(modelSTR.S,2)); zeros(size(modelSTR.S,1), size(modelREF.S,2)) modelSTR.S];
    
    % Set the RHS vector
    modelDelta.b = zeros(size(modelDelta.S,1),1);
    
    % Objective function
    [~,nVars] = size(modelDelta.S);
    % enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    enzymeIds_cindex = enzymeIds + nRxnsREF;
    
    modelDelta.c = zeros(nVars,1);
    modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
    modelDelta.c(enzymeIds) = modelDelta.c(enzymeIds)/etotREF;
    modelDelta.c(enzymeIds_cindex) = -ones(length(modelSTR.rxns(enzymeIds)),1);
    
    % Solve the problem
    LPsolution = solveCobraLP(modelDelta);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full(length(modelREF.rxns)+1:end);
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

    results_CHEM_V.Conditions(k) = growth_conditions(k);
    results_CHEM_V.Pearson(k) = cor{2,1};
    results_CHEM_V.pvalue_P(k) = cor{3,1};
    results_CHEM_V.Spearman(k) = cor{2,2};
    results_CHEM_V.pvalue_S(k) = cor{3,2};
    results_CHEM_V.RMdSE(k) = result;
    results_CHEM_V.NumEs(k) = numEs;

end

%% Loop over growth conditions using GLC_CHEM_mu_0_21_P as reference
% pcGEM integrated with experimental proteomics for reference condition 
load('./Models/GLC_CHEM_mu_0_21_P.mat')                  
modelREF = ecModelP;
clear ecModelP

results_CHEM_P = cell2table(cell(0,0));
% corrVals_CHEM_P = cell2table(cell(0,0));

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(19:26);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('./Models/eciML1515_batch.mat')                    
    modelSTR = ecModel_batch;
    clear ecModel_batch

    % pcGEM integrated with experimental proteomics for suboptimal condition
    modeldir = "./Models/" + current_condition + ".mat";
    load(modeldir)                
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline/Baseline_" + current_condition + ".mat";
    load(baselinedir)  
    baseline = Ecoli;
    clear Ecoli
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(9);

    etotREF = ptotREF * f_REF * sigma_REF;
    etotSTR = ptotSTR * f_STR * sigma_STR;

    modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % Construct the LP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsSTR] = size(modelSTR.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
    [enzymeNum,~] = size(enzymeIds);
    
    % infNum = sum(isinf(modelREF.ub(enzymeIds)));
    % infBounds = find(isinf(modelREF.ub));
    % infBoundsProt = intersect(infBounds, enzymeIds);
    % 
    % etotShare = etotSTR/infNum;
    % modelREF.ub(infBoundsProt) = etotShare;
    
    modelDelta = struct();
    
    % Set the lower and upper bounds for the fluxes
    modelDelta.lb = [modelREF.lb; modelSTR.lb];
    modelDelta.lb(modelDelta.lb==-Inf) = -1000;
    
    modelDelta.ub = [modelREF.ub; modelSTR.ub];
    modelDelta.ub(modelDelta.ub==Inf) = 1000;
    
    % Concatenate the stoichiometric matrices and S-matrices of the two models
    modelDelta.S = [modelREF.S zeros(size(modelREF.S,1), size(modelSTR.S,2)); zeros(size(modelSTR.S,1), size(modelREF.S,2)) modelSTR.S];
    
    % Set the RHS vector
    modelDelta.b = zeros(size(modelDelta.S,1),1);
    
    % Objective function
    [~,nVars] = size(modelDelta.S);
    % enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    enzymeIds_cindex = enzymeIds + nRxnsREF;
    
    modelDelta.c = zeros(nVars,1);
    modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
    modelDelta.c(enzymeIds) = modelDelta.c(enzymeIds)/etotREF;
    modelDelta.c(enzymeIds_cindex) = -ones(length(modelSTR.rxns(enzymeIds)),1);
    
    % Solve the problem
    LPsolution = solveCobraLP(modelDelta);
    
    % Get the solution(s)
    LPsolution.x = LPsolution.full(length(modelREF.rxns)+1:end);
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

    results_CHEM_P.Conditions(k) = growth_conditions(k);
    results_CHEM_P.Pearson(k) = cor{2,1};
    results_CHEM_P.pvalue_P(k) = cor{3,1};
    results_CHEM_P.Spearman(k) = cor{2,2};
    results_CHEM_P.pvalue_S(k) = cor{3,2};
    results_CHEM_P.RMdSE(k) = result;
    results_CHEM_P.NumEs(k) = numEs;

end

%% Merge all results tables
results_table = [results_BATCH; 
    results_CHEM_S; 
    results_CHEM_V; 
    results_CHEM_P];

%%
toc