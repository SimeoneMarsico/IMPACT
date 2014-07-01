function [moduleSet_out,moduleProfile_out,moduleoligosMemb_out]=evalNeighModule(AllNeigh,CurrentPattern,CurrentOligoSymb,GeneProfiles,GeneOligoSymb,thrCorr,num_NeighProf,modSize,simil_type)

moduleSet_out=[];
moduleProfile_out=[];
moduleoligosMemb_out=[];
for j=1:length(AllNeigh)
    
    globalMatrix=[CurrentPattern; GeneProfiles{AllNeigh(j)}]; %matrix numprof x 25
    if(simil_type==0)
        globalCorr=corr(globalMatrix');
    else
        globalCorr=corr_dist(globalMatrix',0);
    end
    %globalOligoSymb=[CurrentOligoSymb;GeneOligoSymb{AllNeigh(j)}];
    globalOligoSymb=[GeneOligoSymb{AllNeigh(j)}];
    upperGlobalCorr=triu(globalCorr,1);
%     correlatedPROF_up=find(upperGlobalCorr(1,:)> thrCorr);
%     correlatedPROF_down=find(upperGlobalCorr(1,:)<-(thrCorr-0.05));
    correlatedPROF_up=find(upperGlobalCorr(1,:)> thrCorr);
    correlatedPROF_down=find(upperGlobalCorr(1,:)<-thrCorr);
    
    
    % updating of current Module
    if (num_NeighProf >= 1)
        len_up = length(correlatedPROF_up);
        len_down = length(correlatedPROF_down);
    else
        len_up = length(correlatedPROF_up) / (size(upperGlobalCorr,2)-1);
        len_down = length(correlatedPROF_down) / (size(upperGlobalCorr,2)-1);
    end
    if len_up > 0 & len_down > 0
        
        if len_up > len_down & len_up >= num_NeighProf
            profile=globalMatrix(correlatedPROF_up,:);
            oligoSymb=globalOligoSymb(correlatedPROF_up-1);
            
            moduleSet_out=[moduleSet_out AllNeigh(j)];
            moduleProfile_out=[moduleProfile_out;profile];
            moduleoligosMemb_out=[moduleoligosMemb_out;oligoSymb];
            clear profile oligoSymb
        elseif len_up < len_down & len_down >= num_NeighProf
            profile=globalMatrix(correlatedPROF_down,:);
            oligoSymb=globalOligoSymb(correlatedPROF_down-1);
            
            moduleSet_out=[moduleSet_out AllNeigh(j)];
            moduleProfile_out=[moduleProfile_out;profile];
            moduleoligosMemb_out=[moduleoligosMemb_out;oligoSymb];
            clear profile oligoSymb
        else
            moduleSet_out=moduleSet_out;
            moduleProfile_out=moduleProfile_out;
            moduleoligosMemb_out=moduleoligosMemb_out;
        
        end
        
    elseif isempty(correlatedPROF_down) & len_up >= num_NeighProf
        profile=globalMatrix(correlatedPROF_up,:);
        oligoSymb=globalOligoSymb(correlatedPROF_up-1);
        moduleSet_out=[moduleSet_out AllNeigh(j)];
        moduleProfile_out=[moduleProfile_out;profile];
        moduleoligosMemb_out=[moduleoligosMemb_out;oligoSymb];
        clear profile oligoSymb
            
    elseif isempty(correlatedPROF_up) & len_down >= num_NeighProf
        profile=globalMatrix(correlatedPROF_down,:);
        oligoSymb=globalOligoSymb(correlatedPROF_down-1);
        moduleSet_out=[moduleSet_out AllNeigh(j)];
        moduleProfile_out=[moduleProfile_out;profile];
        moduleoligosMemb_out=[moduleoligosMemb_out;oligoSymb];
        clear profile oligoSymb
        
            
    else
        moduleSet_out=moduleSet_out;
        
        moduleProfile_out=moduleProfile_out;
        moduleoligosMemb_out=moduleoligosMemb_out;
    end
    
    if length(moduleSet_out)+modSize>100
        break
    end
        
    clear correlatedPROF globalMatrix globalCorr upperGlobalCorr
end
moduleoligosMemb_out;




