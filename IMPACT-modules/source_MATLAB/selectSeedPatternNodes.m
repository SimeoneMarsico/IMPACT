function [IndexSimilarProfiles]=selectSeedPatternNodes(SimilarityMatrix,IndicesProfiles,corrThrSeed)

SeedNodeWindow=find(IndicesProfiles==1);

for i=1:length(SeedNodeWindow)
    
    SimMatrixInd=SeedNodeWindow(i);
    
    [CorrValuesProfiles tempNodesInd]=find(abs(SimilarityMatrix(SimMatrixInd,:))>corrThrSeed);
    
    if isempty(tempNodesInd)
        
        IndicesProfilesCorrelated{i}=[];
        IndicesNodesCorrelated{i}=[];

        signSimilatiries{i}=[];
        valuesSimilatiries{i}=[];
        MembershipProfiles{i}=[];
        
        NumCorrelatedProf(i)=0;        
        CorrSum(i)=0;
        
        continue
    end
    
    IndicesProfilesCorrelated{i}=tempNodesInd;
    IndicesNodesCorrelated{i}=IndicesProfiles(tempNodesInd);
    
    NumCorrelatedProf(i)=length(tempNodesInd);        
    
    signSimilatiries{i}=sign(SimilarityMatrix(SimMatrixInd,tempNodesInd));
    valuesSimilatiries{i}=SimilarityMatrix(SimMatrixInd,tempNodesInd);
    
    MembershipProfiles{i}=IndicesProfiles(tempNodesInd);
    
    CorrSum(i)=sum(abs(SimilarityMatrix(SimMatrixInd,tempNodesInd)));
        
    clear CorrValuesProfiles tempNodesInd
end

clear k

for i=1:length(SeedNodeWindow)
    IndSeed=SeedNodeWindow(i);
    
    if isempty(IndicesProfilesCorrelated{IndSeed})
        SimilarProfIndEnd(i)=NumCorrelatedProf(i);
        continue
    end
    %check if profiles from the same node have at the same time positive
    %and negatie sign
    temp_sign=signSimilatiries{IndSeed};
    temp_membership=MembershipProfiles{IndSeed};
    
    unique_membership=unique(temp_membership);
    if length(temp_membership)==length(unique_membership)
        SimilarProfIndEnd=NumCorrelatedProf;
        break
    end
    
    
    for j=1:length(unique_membership)
        indNOde=find(temp_membership==unique_membership(j));
        prodSignNOde=prod(sign(temp_sign(indNOde)));
        
        if prodSignNOde==-1
            flag_sign=0;%this flag allow to check if there is among the correlated profiles, 
                        %some belonging to the same node that have an opposite sign. 
                        %flag_sign==0 means that this is the case; flag_sign==1 means that this case is note present
            continue
        end
        flag_sign=1;
    end
    
    if flag_sign==1
        SimilarProfIndEnd(i)=NumCorrelatedProf(i);
    elseif flag_sign==0
        SimilarProfIndEnd(i)=0;
    end
    
    clear temp_sign temp_membership indNOde sumSignNOde unique_membership
end
if max(SimilarProfIndEnd)==0
    IndexSimilarProfiles=[];
    return

end

SimilarProfInd=find(SimilarProfIndEnd==max(SimilarProfIndEnd));

if isempty(SimilarProfInd)
    IndexSimilarProfiles=[];
elseif length(SimilarProfInd)>1
    [MaxCorr,IndMaxCorr]=max(CorrSum(SimilarProfInd));
%    IndexSimilarProfiles=[SimilarProfIndEnd(SimilarProfInd(IndMaxCorr)) IndicesProfilesCorrelated{SimilarProfInd(IndMaxCorr)}];
    IndexSimilarProfiles=[SimilarProfInd(IndMaxCorr) IndicesProfilesCorrelated{SimilarProfInd(IndMaxCorr)}];
else
%    IndexSimilarProfiles=[SimilarProfIndEnd(SimilarProfInd) IndicesProfilesCorrelated{SimilarProfInd}]; 
    IndexSimilarProfiles=[SimilarProfInd IndicesProfilesCorrelated{SimilarProfInd}];
end
