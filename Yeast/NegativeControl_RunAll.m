%
% Retrieves kcat
%
% Code partially adapted from MOMA.m (COBRA Toolbox 3)
%
%% Cleaning the workspace and the command window
clear;clc
tic
changeCobraSolver('gurobi', 'QP');

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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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

ptotREF = 0.265;

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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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
    modelSTR_Exp = ecModelP;
    clear ecModelP

    % Baseline for calculating comparisons
    baselinedir = "./Baseline2/Baseline_Norm2_" + current_condition + ".mat";
    load(baselinedir) 
    
    % Retrieving kcat values
    proteins = modelSTR.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    kcat = [];
    
    for i=1:numel(proteins)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(proteins(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
        kcat(i,1) = min(-1./modelSTR.S(pIdx(i),rxnIdx{i}));
    end
    
    % Calculate correlation with experimental values
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    
    pred = {};
    pred(:,1) = modelSTR.rxns(enzymeIds);
    pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
    pred(:,2) = num2cell(kcat(:,1));
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
    rmse.Predicted = log10(abs(kcat));
    
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

% corrVals_final = [corrVals_LAHTVEE;
%     corrVals_YU2021_N30;
%     corrVals_YU2021_Gln;
%     corrVals_YU2021_Phe;
%     corrVals_YU2021_Ile;
%     corrVals_YU2020];

%%
toc