function pvalue = pvalue_fromFrequencies(referenceProfile,actualCase,randomizationMatrix)

Nrand=size(randomizationMatrix,1);
for i=1:length(referenceProfile)
    
    if isempty(referenceProfile{i})
        pvalue(i)=1;
        
    else
        pvalue(i)=(length(find(randomizationMatrix(:,i)>=actualCase(i)))+0)/Nrand;
    end
end
