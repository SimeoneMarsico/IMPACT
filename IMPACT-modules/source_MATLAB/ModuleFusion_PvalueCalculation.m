function [outputFlag] = ModuleFusion_PvalueCalculation(ModuleResultsFile,EstimatedPatternProbabilities,FinalResultsFile,SIMIL_THR,N_MAX,simil_type)
% function [outputFlag] = ModuleFusion_PvalueCalculation(moduleResultsFile,EstimatedPatternProbabilities,FinalResultsFile,SIMIL_THR)
% This function requires as inputs:
% 1 - filename of mat file storing the results obtained from the module analysis routine
% 2 - filename of mat file storing the empirically estimated probabilities that are profiles specific (mat file)
% 3 - filename of the mat file that will contain the final results (after
% module fusion step and pvalue estimation)
% 4 - similarity threshold 


load(ModuleResultsFile);

load(EstimatedPatternProbabilities);

P_PseudoSeedProfiles = P_SeedPatternUsed; clear P_SeedPatternUsed
PseudoSeedProfiles = SeedPatternUsed; clear Estimated_patternProbabilities
clear Estimated_patternProbabilities

modulePvalue=gene_network_random(AllPROFnet,AdjMatrixNet_Screen,SeedNodeUsed,SeedPatternUsed,PseudoSeedProfiles,moduleSet,P_PseudoSeedProfiles,SIMIL_THR,N_MAX,simil_type);


%code for swapping the profiles in case that they are anticorrelated
SwappedProf={};

h = waitbar(0,'swapping profiles...');
for i=1:length(moduleProfile)
    waitbar(i/length(moduleProfile));
    if(simil_type==0)
        corrWithSeedProf=corr([SeedPatternUsed(i,:);moduleProfile{i}]');
    else
        corrWithSeedProf=corr_dist([SeedPatternUsed(i,:);moduleProfile{i}]',0);
    end
    signCorr=sign(corrWithSeedProf(1,2:end));
    indProfToSwap=find(signCorr==-1);
    temp=moduleProfile{i}(indProfToSwap,:)*-1;
    SwappedProf{i}=moduleProfile{i};
    SwappedProf{i}(indProfToSwap,:)=temp;
    moduleReferenceProfile(i,:)=mean(SwappedProf{i})';
    
end
if (exist('h'))
    close(h);
end

%code for fuse modules on the basis of the correlation between seed
%pattern and the overlapping between modules
%for checking the correlation criterion on the module profile instead of
%the module seed profile change: corrSeedProf --> corrModProf
% newModule1={};
% newUnionINdices1=[];VectorJ=[];
% k=1;
% h = waitbar(0,'fusing modules...');
% for i=1:length(moduleSet)
%     waitbar(i/length(moduleSet));
% 
%     if ~ismember(i,VectorJ)
%         indicesTemp=i;
%         mod1=moduleSet{i};
% 
%         for j=i+1:length(moduleSet)
% 
%             if ~ismember(j,VectorJ)
% 
%                 mod2 = moduleSet{j};
%                 [commonNodes, indCommonNodes1, indCommonNode2]=intersect(mod1, mod2);
% 
%                 %corrModProf=corr(moduleReferenceProfile(i,:)',moduleReferenceProfile(i,:)');
%                 corrModProf=corr2_dist(moduleReferenceProfile(i,:)',moduleReferenceProfile(j,:)',0);
%                 
%                 if commonNodes
% 
%                     if corrModProf>=SIMIL_THR && (length(commonNodes)/length(mod1))> SIMIL_THR || corrModProf>=SIMIL_THR && (length(commonNodes)/length(mod2))> SIMIL_THR  
% 
%                        if length(mod1)>40 || length(mod2)>40
%                             continue
%                         end
%                         mod1 = union(mod1,mod2);
%                         fMod = mod1;
%                         indicesTemp = [indicesTemp j];
%                         VectorJ = [VectorJ j];
%                     end
%                 end
%         
%             end
%         
%         end
%         if exist('fMod')
%             newModule1{k}=fMod;
%             newUnionINdices1{k}=indicesTemp;
%             k=k+1;            
%             clear indicesTemp fMod
%         end
%     end
% end
% if (exist('h'))
%     close(h);
% end

%at the end of the checking phase, newModule1 contains the new modules
%obtainend after merging the pre-existing modules found with the method
%and neUnionINdices1 contains the id of the original modules used for the
%mergin.
%second step union _ consensus modules detection! == necessary for the
%pvalue computation!
% NewUnionIndices_vector=cell2mat(newUnionINdices1);
% RemainingModule=setdiff(1:length(moduleProfile),NewUnionIndices_vector);
% 
% clear i j
% 
% %newModules after fusion
% for i=1:length(newUnionINdices1)
%     
%     if length(newUnionINdices1{i})>2
%         temp=[];
%         for j=1:length(newUnionINdices1{i})
%             temp=union(temp,moduleSet{newUnionINdices1{i}(j)});
%         end
%         newModulesSet_first{i}=temp;
%     elseif length(newUnionINdices1{i})<=2
%         newModulesSet_first{i}=union(moduleSet{newUnionINdices1{i}});
%     end
%     clear temp
%            
% end
% 
% for j=1:length(RemainingModule)
%     newModulesSet_second{j}=moduleSet{RemainingModule(j)};
%     SeedNodeUsedNew(j)=SeedNodeUsed(RemainingModule(j));
%     SeedPatternUsedNew(j,:)=SeedPatternUsed(RemainingModule(j),:);
% end
% 
% %create new variables for computing the new pvalues
% SeedNodeUsedNew_second=[];
% SeedPatternUsedNew_second=[];counter=1;moduleSetNEW_forPValue={};
% for i=1:length(newUnionINdices1)
%     for j=1:length(newUnionINdices1{i})
%         
%         SeedNodeUsedNew_second(counter)=SeedNodeUsed(newUnionINdices1{i}(j));
%         SeedPatternUsedNew_second(counter,:)=SeedPatternUsed(newUnionINdices1{i}(j),:);
%         moduleSetNEW_forPValue{counter}=newModulesSet_first{i};
%         counter=counter+1;
%     end
% end
% 
% NEWModuleSet=[moduleSetNEW_forPValue newModulesSet_second];
% NEWSeedNodeUsed=[SeedNodeUsedNew_second SeedNodeUsedNew];
% NEWSeedPatternUsed=[SeedPatternUsedNew_second;SeedPatternUsedNew];
% 
% %% MODULE PVALUE CALCULATION 
% modulePvalueNEW=gene_network_random(AllPROFnet,AdjMatrixNet_Screen,NEWSeedNodeUsed,NEWSeedPatternUsed,PseudoSeedProfiles,NEWModuleSet,P_PseudoSeedProfiles,SIMIL_THR,N_MAX);

% VectorPvalues2=[];
% VectorPvalues1=[];
% FinalFusedModulePvalue=[];
% sizeModules1=[];
% sizeModules2=[];
% overlapOligoid=[];
% indicesOriginal=[];
% indicesNEW=[];
% counter=1;
% h = waitbar(0,'pvalue calculation...');
% for i = 1:length(newUnionINdices1)
%     waitbar(i/length(newUnionINdices1));
%     temp = modulePvalue(newUnionINdices1{i});
%     VectorPvalues1 = [VectorPvalues1 temp];
%    
%     l1 = length(newUnionINdices1{i});
%     tempOligoId = moduleOligoSymb{newUnionINdices1{i}(1)};
%     for j=1:l1
%         
%         indTemp1 = find(NEWSeedNodeUsed==SeedNodeUsed(newUnionINdices1{i}(j)));
%         temp2(j) = modulePvalueNEW(indTemp1);
%         temp3 = moduleOligoSymb{newUnionINdices1{i}(j)};
%         sizeModules1 = [sizeModules1 length(moduleSet{newUnionINdices1{i}(j)})];
%         indicesOriginal = [indicesOriginal newUnionINdices1{i}(j)];
%         indicesNEW = [indicesNEW i];
%         overlapOligoid = [overlapOligoid (length(intersect(tempOligoId,temp3))/min(length(tempOligoId),length(temp3)))];
%         tempOligoId = union(tempOligoId,temp3);
%         clear indTemp1
%     end
%     sizeModules2 = [sizeModules2 repmat(length(newModule1{i}),1,length(newUnionINdices1{i}))];
%     VectorPvalues2 = [VectorPvalues2 temp2];
%     
%     FinalFusedModulePvalue = [FinalFusedModulePvalue min(temp2)];
%     
%     clear temp temp2 l1 temp3 tempOligoId
% end
% if (exist('h'))
%     close(h);
% end
% %create the final list of pvalues
% %criteria:
% % 1) if a module is the complete subset of another bigger one, the pvalue
% % of this module is the pvalue of the biggest
% % 2) if 2 modules are aprtially overlapping, the pvalue is the mean among
% % all.
% % 3) the min pvalue obtained testing the different seed.
% 
% moduleProfileFused={};
% moduleProfileFusedtemp=[];
% moduleSetFused={};
% moduleSetFusedtemp=[];
% moduleOligoSymbFused={};
% moduleOligoSymbFusedtemp=[];
% SeedNodeUsedFused=[];
% SeedPatternUsedFused=[];
% modulePvalueFused=[];
% for i=1:length(newUnionINdices1)
%     l1=length(newUnionINdices1{i});
%     
%     for j=1:l1
%         
%         indTemp1(j)=find(NEWSeedNodeUsed==SeedNodeUsed(newUnionINdices1{i}(j)));
%         temp2(j)=modulePvalueNEW(indTemp1(j));
%         
%         moduleProfileFusedtemp=[moduleProfileFusedtemp;moduleProfile{newUnionINdices1{i}(j)}];
%         moduleSetFusedtemp=[moduleSetFusedtemp moduleSet{newUnionINdices1{i}(j)}];
%         moduleOligoSymbFusedtemp=[moduleOligoSymbFusedtemp;moduleOligoSymb{newUnionINdices1{i}(j)}];
%         
%     end
%     
%     moduleProfileFusedtemp2=[];
%     [moduleOligoSymbFusedtemp2 indUniqueOligo]=unique(moduleOligoSymbFusedtemp);
%     moduleProfileFusedtemp2=moduleProfileFusedtemp(indUniqueOligo,:);
%     
%     [minPval, indminPval]=min(temp2);
%     SeedNodeUsedFused(i)=NEWSeedNodeUsed(indTemp1(indminPval));
%     %SeedPatternUsedFused(i,:)=NEWSeedPatternUsed(indTemp1(indminPval),:);
%     SeedPatternUsedFused(i,:)=median(SeedPatternUsed(newUnionINdices1{i},:));
%     moduleProfileFused{i}=moduleProfileFusedtemp2;
%     moduleOligoSymbFused{i}=moduleOligoSymbFusedtemp2;
%     moduleSetFused{i}=unique(moduleSetFusedtemp);
%     modulePvalueFused(i)=minPval;
%     clear temp2 indTemp minPval indminPval moduleOligoSymbFusedtemp2 moduleProfileFusedtemp2
%     moduleProfileFusedtemp=[]; moduleOligoSymbFusedtemp=[]; moduleSetFusedtemp=[];
% end
% clear moduleOligoSymbFusedtemp moduleProfileFusedtemp moduleSetFusedtemp
% %append old Results (all those modules that have not been fused)
% 
% k=length(moduleProfileFused);
% for i=1:length(RemainingModule)
%     indRem=RemainingModule(i);
%     SeedNodeUsedFused(k)=SeedNodeUsed(indRem);
%     SeedPatternUsedFused(k,:)=SeedPatternUsed(indRem,:);
%     moduleProfileFused{k}=moduleProfile{indRem};
%     moduleOligoSymbFused{k}=moduleOligoSymb{indRem};
%     moduleSetFused{k}=moduleSet{indRem};
%     modulePvalueFused(k)=modulePvalue(indRem);
%     k=k+1;clear indRem
% end


%% SAVE RESULTS INTO A RESULTS FILE

save(FinalResultsFile,'modulePvalue*','moduleProfile*','moduleSet*','moduleOligoSymb*','SeedNodeUsed*','SeedPatternUsed*','modulePvalue*','AdjMatrixNet_Screen','AllGenesNet_Screen','AllPROFnet', 'SwappedProf', 'moduleReferenceProfile');

if isempty(moduleProfile)
    outputFlag = 0;
else
    outputFlag=1;
end