function writeProtCounts(condition,prot_parameters,name)
initial_prots  = length(prot_parameters{1});
filtered_prots = length(prot_parameters{2});
matched_prots  = length(prot_parameters{3});
model_prots    = length(prot_parameters{4});
mass_coverage  = prot_parameters{5};
disp(['The total number of proteins in the dataset is:                ' num2str(initial_prots)])
disp(['The total number of proteins in the filtered dataset is:       ' num2str(filtered_prots)])
disp(['The total number of filtered proteins present in the model is: ' num2str(matched_prots)])
disp(['The mass ratio between measured and unmeasured protein is:     ' num2str(mass_coverage)])

T = table(initial_prots,filtered_prots,matched_prots,model_prots,mass_coverage);
writetable(T,['../../models/prot_constrained/' name '/prot_counts_' condition '.txt'],'Delimiter','\t')
end