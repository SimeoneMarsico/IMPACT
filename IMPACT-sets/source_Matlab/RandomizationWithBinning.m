function [Rand_selectedProfiles Rand_SelectedGenes]=RandomizationWithBinning(newComplexes,newBinningSetProt,Nrandomizations,vectRandTemp,thresholdSimil,N, Realselected,simil_type)
% for randomizations
ProvaBinning_GlobalBestS_1000 = zeros(Nrandomizations,length(newComplexes));
ProvaBinning_GlobalRandomSelGenes_10000 = zeros(Nrandomizations,length(newComplexes));
h = waitbar(0,'running randomizations...');
for i=1:length(newComplexes)
    i
    if (isempty(Realselected{i}))
        continue;
    end
    waitbar(i/length(newComplexes));
    if ~isempty(newComplexes{i})
        currComplex=newComplexes{i};
        sizeComplex=length(currComplex);
        bestS=zeros(Nrandomizations,1);
        for j=1:Nrandomizations

            for k=1:sizeComplex
                currProt=currComplex{k};        % selection of the protein in the current complex
                numProfProt=size(currProt,1);   % read the number of oligo profiles available for the protein
                temp_ind = find(numProfProt <= N);
                if(~isempty(temp_ind))
                    numProfProt_choice = temp_ind(1);
                else
                    numProfProt_choice = length(N);
                end
                if (sizeComplex == 1)
                    a = newBinningSetProt{numProfProt_choice};
                    b = cell2mat(a);
                    randomIndexProtein=randsample(size(b,1),numProfProt);
                    RandM{k}=b(randomIndexProtein,:);
                else
                    randomIndexProtein=randsample(vectRandTemp(numProfProt_choice,1),1);
                    RandM{k}=newBinningSetProt{numProfProt_choice}{randomIndexProtein};
                end
                clear currProt numProfProt randomIndexProtein
            end

            randComplex=RandM;
            [bestSize,numGenesSelected]=searchigEnrichedPRofiles2Noprofiles(randComplex,thresholdSimil,simil_type);
            bestS(j)=bestSize;
            RandSelgenes(j)=numGenesSelected;
            clear selectedProfiles referenceProfile bestSize randprot numGenesSelected
            clear randComplex RandM k 

        end
        clear currComplex sizeComplex
        ProvaBinning_GlobalBestS_1000(:,i) = bestS'; 
        ProvaBinning_GlobalRandomSelGenes2_1000(:,i) = RandSelgenes'; 
        clear bestS RandSelgenes RandSelgenes    
    end
end
if(exist('h'))
    close(h);
end
Rand_selectedProfiles=ProvaBinning_GlobalBestS_1000;
Rand_SelectedGenes=ProvaBinning_GlobalRandomSelGenes2_1000;

     