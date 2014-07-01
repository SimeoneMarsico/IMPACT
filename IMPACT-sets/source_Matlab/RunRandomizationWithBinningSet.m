function [outputFlag]=RunRandomizationWithBinningSet(inputFile,outputFile,numRandomization,ThresholdRand, N,Realselected,simil_type)
% function [ouputFlag]=RunRandomizationWithBinning(inputFile,outputFile) 
% inputFile is a mat file containing all necessary information
% outputFlag = procedure finished

if ischar(numRandomization)
   numRandomization=str2num(numRandomization);
end

if ischar(ThresholdRand)
   ThresholdRand=str2double(ThresholdRand);
end

load(inputFile)
% input file contains: SetProfiles Bins vectRandTemp 
[Rand_selectedProfiles Rand_SelectedGenes]=RandomizationWithBinning(SetProfiles,Bins,numRandomization,vectRandTemp,ThresholdRand,N,Realselected,simil_type);
save(outputFile,'Rand_selectedProfiles','Rand_SelectedGenes');
outputFlag=1;
