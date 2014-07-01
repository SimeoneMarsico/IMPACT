function [moduleSet1,moduleProfile1,moduleOligoSymb1,moduleScore1,SeedNodeUsed1,SeedPatternUsed1]=TOPModule(PotentialSeedNodeVect,SeedNodePatternMatrix,SeedNodeSymbols,SeedNodeProfilesSelected,gene_profile1,oligoSymbols,AdjMatrix,CorrT,typeModScore,minNumProfiles,simil_type)

fprintf('\n Searching procedure started ... \n')
[moduleSet1,moduleProfile1,moduleOligoSymb1,moduleScore1,SeedNodeUsed1,SeedPatternUsed1]=pipelineSearchingProcedure2(PotentialSeedNodeVect,SeedNodePatternMatrix,SeedNodeSymbols,SeedNodeProfilesSelected,AdjMatrix,gene_profile1,oligoSymbols,CorrT,typeModScore,minNumProfiles,simil_type);
fprintf('\n\t\t ... searching procedure completed successfully! \n\n')



