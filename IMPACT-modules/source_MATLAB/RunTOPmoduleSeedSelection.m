function [PotentialSeedNodeDef,SeedNodePatternOriginal,SeedNodeProfilesOriginal,SeedNodeOligoSymbols,SeedNodeProfilesOriginal2,SeedNodeOligoSymbols2]=RunTOPmoduleSeedSelection(networkFile,ExperimentalDataFile,typeSeedNodeSelection,seedPatternThr,OutputResultFile,minNumProfiles,simil_type)
% networkFile: mat file containing the adjacency matrix

%load info on network and experimental data
load(networkFile); % it contains : AdjMatrixNet and AllGenesNet
%
load(ExperimentalDataFile)
%allOligoEns
%allOligoProf
%allOLIGOName

%control on the numeric parameters

if ischar(seedPatternThr)
    seedPatternThr=str2num(seedPatternThr);
end

if ischar(typeSeedNodeSelection)
    typeSeedNodeSelection=str2num(typeSeedNodeSelection);
end

if ischar(minNumProfiles)
    minNumProfiles=str2num(minNumProfiles);
end

%mapping experimental data onto the network
[AllGenesNet_Screen,indNet,indScreen]=intersect(AllGenesNet,allOligoEns);
AdjMatrixNet_Screen=AdjMatrixNet(indNet,indNet);
AllPROFnet=allOligoProf(indScreen);
AllOligoSymb=allOLIGOName(indScreen);

for i=1:length(AllPROFnet)
	numProfilesAllNet(i)=size(AllPROFnet{i},1);
end
% computation of the seed nodes
ADJ_size=length(AdjMatrixNet_Screen);

% instruction in case of oligo profiles
%[PotentialSeedNodeOriginal]=SelectSeedNodes(MedianCorrelation1,Nseed,Nprof,num_profile,ADJ_size);
%instruction in case of gene profiles
NumProfMatrix=numProfilesAllNet(1:ADJ_size);
PotentialSeedNodeOriginal=find(NumProfMatrix>=minNumProfiles);
[PotentialSeedNodeDef,SeedNodePatternOriginal,SeedNodeProfilesOriginal,SeedNodeOligoSymbols,SeedNodeProfilesOriginal2,SeedNodeOligoSymbols2]=SeedsSelectionFunc(AdjMatrixNet_Screen,PotentialSeedNodeOriginal,AllPROFnet,AllOligoSymb,typeSeedNodeSelection,seedPatternThr,NumProfMatrix,simil_type);

save(OutputResultFile, 'PotentialSeedNodeDef','SeedNodePatternOriginal','SeedNodeProfilesOriginal','SeedNodeOligoSymbols','SeedNodeProfilesOriginal2','SeedNodeOligoSymbols2');