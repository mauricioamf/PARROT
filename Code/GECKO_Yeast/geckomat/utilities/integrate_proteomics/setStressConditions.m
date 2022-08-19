function condModel = setStressConditions(model,Drate,pos,expData,NGAM,interval)
minProt = true;
ecFlag  = true;
%Rescale biomass set fitted GAM for std conditions and fit non-growth
%associated maintenance for stress conditions
cd ..
condModel = setChemostatConstraints(model,pos,Drate,minProt,0.01);
cd integrate_proteomics
disp(' ')
condModel = fitNGAM(condModel,NGAM,expData,interval,ecFlag);
end