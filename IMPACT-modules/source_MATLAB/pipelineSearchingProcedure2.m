function [moduleSet_all,moduleProfile_all,moduleOligoSymbol_all,moduleScore_all,SeedNodeUsed,SeedPatternUsed]=pipelineSearchingProcedure2(SeedNodeList,SeedNodePattern,SeedNodeSymbols,SeedNodeProfilesSelected,AM,GeneProfiles,GeneOligoSymb,CorrThr,typeModScore,num_NeighProf,simil_type)

moduleSet_all={};
moduleProfile_all={};
moduleOligoSymbol_all={};
SeedNodeUsed=[];
moduleScore_all=[];
SeedPatternUsed=[];
ii=1;

h = waitbar(0,'expanding seed nodes');
for i=1:length(SeedNodeList)
    waitbar(i/length(SeedNodeList));
    
    CurrentPattern=SeedNodePattern(i,:); % 1 x 25  %selection of the module pattern
    moduleProfile=SeedNodeProfilesSelected{i};%CurrentPattern;
    moduleSet=[SeedNodeList(i)]; %initialization of module: the module contains the seed node
    moduleSymbol=[SeedNodeSymbols{i}];
    CurrentPatternOligoSymb=moduleSymbol;
    flag_recursion=1;

    while flag_recursion == 1 % flag that determines the end of the for loop
        AllNeigh=searchNeighModule_improved(AM,moduleSet);
        AllNeigh=setdiff(AllNeigh,moduleSet);
        if isempty(AllNeigh)
            break
        end
        thrCorr=CorrThr;
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        moduleSize=length(moduleSet);
        [neighInd,neighProfile,neighOligoSymb]=evalNeighModule(AllNeigh,CurrentPattern,CurrentPatternOligoSymb,GeneProfiles,GeneOligoSymb,thrCorr,num_NeighProf,moduleSize, simil_type);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(moduleSet)>50 % max module size 
            break
        end
        if isempty(neighInd)
            flag_recursion=0;
            continue
        else

            moduleSet=[moduleSet neighInd];
            moduleProfile=[moduleProfile; neighProfile];
            moduleSymbol=[moduleSymbol;neighOligoSymb];
            flagrecursion=1;
        end
        
        
    end
%   UPDATE OF OUTPUT VARIABLES
    
    SeedNodeUsed(ii)=SeedNodeList(i);
    SeedPatternUsed(ii,:)=CurrentPattern;
    moduleSet_all{ii}=moduleSet;
    moduleProfile_all{ii}=moduleProfile;
    moduleOligoSymbol_all{ii}=moduleSymbol;
    if(simil_type==0)
        score_temp=triu(corr(moduleProfile'),1);
    else
        score_temp=triu(corr_dist(moduleProfile',0),1);
    end
    moduleScore_all(ii,:)=MedianValueFromMatrix(score_temp);
    ii=ii+1;

    NtempProf=size(moduleProfile,1);
    message1=strcat(' - - Seed node: ', num2str(i), '   module size:', num2str(length(moduleSet)),'   N profiles:', num2str(NtempProf),'\n');
    
    clear moduleSet moduleProfile currentSeedNode AllNeigh message1 score_temp neighOligoSymb
    
end
if (exist('h'))
    close(h);
end
length(SeedNodeList)