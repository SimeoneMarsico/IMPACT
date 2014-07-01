function [PotentialSeedNodeDef,SeedNodePatternOriginal,SeedNodeProfilesOriginal,SeedNodeOligoSymbols,SeedNodeProfilesOriginal2,SeedNodeOligoSymbols2]=SeedsSelectionFunc(AdjMatrix,PotentialSeedNodeOriginal,gene_profile1,geneOligoSymbols,seedPatternSelection,corrThrSeed,num_profile,simil_type)


% seedPatternSelection: can be 1 for having the seed set selection and 0
% for having the seed node selection.

% selection of median pattern
k=1;
PotentialSeedNodeDef=[];
SeedNodePatternOriginal=[];
SeedNodeOligoSymbols=[];

for i=1:length(PotentialSeedNodeOriginal)
    
    if PotentialSeedNodeOriginal(i)==3
        PotentialSeedNodeOriginal;
    end
    if seedPatternSelection==1
        
        NeighSeed = searchNeighModule_improved(AdjMatrix,PotentialSeedNodeOriginal(i));
        
        if isempty(NeighSeed)
            continue
        end

        MatrixProfiles=gene_profile1{PotentialSeedNodeOriginal(i)};
        MatrixOligoSymbols=geneOligoSymbols{PotentialSeedNodeOriginal(i)};
        IndexProfilesmembership=repmat(1,size(gene_profile1{PotentialSeedNodeOriginal(i)},1),1);
        for j=1:length(NeighSeed)
            
            MatrixProfiles=[MatrixProfiles;gene_profile1{NeighSeed(j)}];
            MatrixOligoSymbols=[MatrixOligoSymbols;geneOligoSymbols{NeighSeed(j)}];
            IndexProfilesmembership=[IndexProfilesmembership;repmat(j+1,size(gene_profile1{NeighSeed(j)},1),1)];
            
        end
        if(simil_type==0)
            CorrMatrixProfile=corrcoef(MatrixProfiles');
        else
            CorrMatrixProfile=corrcoef_dist(MatrixProfiles',0);
        end
        
        clear j
        
        for j=1:length(CorrMatrixProfile)
            CorrMatrixProfile(j,j)=0;
        end
%         i
        [IndexSimilarProfilesTest]=selectSeedPatternNodes(CorrMatrixProfile,IndexProfilesmembership,corrThrSeed);
        IndexSimilarProfilesTest=unique(IndexSimilarProfilesTest);
        SeedNodeOligoSymbols_temp=MatrixOligoSymbols(IndexSimilarProfilesTest);
        Seedsymbol= strtok(geneOligoSymbols{PotentialSeedNodeOriginal(i)}(1),':');
        indicesSeedOnlyProfiles=strmatch(Seedsymbol,SeedNodeOligoSymbols_temp);
                
        if ~isempty(IndexSimilarProfilesTest) && ~isempty(indicesSeedOnlyProfiles)
            
            PotentialSeedNodeDef(k)=PotentialSeedNodeOriginal(i);

            SelectedSeedProfiles=MatrixProfiles(IndexSimilarProfilesTest,:); 

            SeedNodePatternOriginal(k,:)=median(SelectedSeedProfiles);
            
            SeedNodeProfilesOriginal{k}=SelectedSeedProfiles;
            
            SeedNodeOligoSymbols{k}=MatrixOligoSymbols(IndexSimilarProfilesTest);
            
            SeedNodeProfilesOriginal2{k}=SelectedSeedProfiles(indicesSeedOnlyProfiles,:);
            SeedNodeOligoSymbols2{k}=SeedNodeOligoSymbols{k}(indicesSeedOnlyProfiles);            
            k=k+1;
            clear SelectedSeedProfiles Seedsymbol indicesSeedOnlyProfiles
        end
        
        
        clear NeighSeed MatrixProfiles MatrixOligoSymbols IndexProfilesmembership CorrMatrixProfile j

    elseif seedPatternSelection==0
        
        PotentialSeedNodeDef(i)=PotentialSeedNodeOriginal(i);
        
        if num_profile(PotentialSeedNodeOriginal(i)) >1

            SeedNodePatternOriginal(i,:)=median(gene_profile1{PotentialSeedNodeOriginal(i)});
            SeedNodeProfilesOriginal{i}=gene_profile1{PotentialSeedNodeOriginal(i)};
        else

            SeedNodePatternOriginal(i,:)=gene_profile1{PotentialSeedNodeOriginal(i)};
            SeedNodeProfilesOriginal{i}=gene_profile1{PotentialSeedNodeOriginal(i)};
        end
        SeedNodeOligoSymbols{i}=geneOligoSymbols{PotentialSeedNodeOriginal(i)};
    end
end