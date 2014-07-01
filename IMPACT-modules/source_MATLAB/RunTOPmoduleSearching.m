function [outputFlag]=RunTOPmoduleSearching(networkFile,ExperimentalDataFile,SeedNodeFile,VectorSeedNodeStartingFile,moduleExpansionThr,OutputResultFile,minNumProfiles,simil_type)
% function []=RunTOPmoduleSearching(networkFile,NetworkGeneListFile,ExperimentalDataFile)
% networkFile: mat file containing the adjacency matrix

%load info on network and experimental data
load(networkFile);
%
load(ExperimentalDataFile)
%allOligoEns
%allOligoProf
%allOLIGOName

load(SeedNodeFile)
DesideredSeedNOdes=dlmread(VectorSeedNodeStartingFile);
%control on the numeric parameters

if length(DesideredSeedNOdes)~=length(PotentialSeedNodeDef)
    fprintf('ERROR --- Seed node vector incopatible with seed nodes loaded from file ---');
    error('Wrong number of vector input')
else
    fprintf('Seed node vector copatible with seed nodes loaded from file');
end


if ischar(moduleExpansionThr)
    moduleExpansionThr=str2num(moduleExpansionThr);
end


if ischar(minNumProfiles)
    minNumProfiles=str2num(minNumProfiles);
end



%mapping experimental data onto the network
[AllGenesNet_Screen,indNet,indScreen]=intersect(AllGenesNet,allOligoEns);
AdjMatrixNet_Screen=AdjMatrixNet(indNet,indNet);
AllPROFnet=allOligoProf(indScreen);
AllOLIGOnet=allOLIGOName(indScreen);
for i=1:length(AllPROFnet)
	numProfilesAllNet(i)=size(AllPROFnet{i},1);
end
% computation of the seed nodes
ADJ_size=length(AdjMatrixNet_Screen);

IndStartingNodes=find(DesideredSeedNOdes==1);
PotentialSeedNodeDef_touse=PotentialSeedNodeDef(IndStartingNodes);
SeedNodePatternOriginal_touse=SeedNodePatternOriginal(IndStartingNodes,:);
SeedNodeSybols_touse=SeedNodeOligoSymbols2(IndStartingNodes);
SeedNodeProfiles_touse=SeedNodeProfilesOriginal2(IndStartingNodes);
[moduleSet,moduleProfile,moduleOligoSymb,moduleScore,SeedNodeUsed,SeedPatternUsed]=TOPModule(PotentialSeedNodeDef_touse,SeedNodePatternOriginal_touse,SeedNodeSybols_touse,SeedNodeProfiles_touse,AllPROFnet,AllOLIGOnet,AdjMatrixNet_Screen,moduleExpansionThr,'r_square',minNumProfiles,simil_type);
if ~isempty(moduleScore)
    outputFlag=1;
else
    outputFlag=0;
end
for i=1:length(moduleSet)
    ModuleSize(i)=length(moduleSet{i});
end

indNotEmtpyMod=find(ModuleSize>1);

moduleSetFinal=moduleSet(indNotEmtpyMod);
moduleProfileFinal=moduleProfile(indNotEmtpyMod);
moduleOligoSymbFinal=moduleOligoSymb(indNotEmtpyMod);
moduleScoreFinal=moduleScore(indNotEmtpyMod,:);
SeedNodeUsedFinal=SeedNodeUsed(indNotEmtpyMod);
SeedPatternUsedFinal=SeedPatternUsed(indNotEmtpyMod,:);
clear moduleSet moduleProfile moduleOligoSymb moduleScore SeedNodeUsed SeedPatternUsed

moduleSet=moduleSetFinal;
moduleProfile=moduleProfileFinal;
moduleOligoSymb=moduleOligoSymbFinal;
moduleScore=moduleScoreFinal;
SeedNodeUsed=SeedNodeUsedFinal;
SeedPatternUsed=SeedPatternUsedFinal;
clear moduleSetFinal moduleProfileFinal moduleOligoSymbFinal moduleScoreFinal SeedNodeUsedFinal SeedPatternUsedFinal

save(OutputResultFile, 'moduleSet','moduleProfile','moduleOligoSymb','moduleScore','SeedNodeUsed','SeedPatternUsed','AdjMatrixNet_Screen','AllGenesNet_Screen','AllPROFnet');