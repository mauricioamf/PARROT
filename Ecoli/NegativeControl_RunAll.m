%
% Retrieves Kcat
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