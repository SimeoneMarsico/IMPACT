function IMPACT_modules(do_params, alg_params)
% INPUT:
% do_params = [do_seed_nodes do_module_search do_randomizations do_statistics]
% where do_seed_nodes     = binary, [0,1], do the import (convert input text files into matlab environment files)
%       do_module_search  = binary, [0,1], do the search procedure on sets
%       do_randomizations = binary, [0,1], do the randomizations
%       do_statistics     = binary, [0,1], calculate p-values and export results based on the randomization results
% alg_params = [SIMIL_THR MIN_NUMB_PROF N_RAND N_MAX]
% where SIMIL_THR = number, [0-1], similarity threshold (Pearson correlation coefficient) for set search and randomizations
%       MIN_NUMB_PROF = number, if >= 1 it is the number of profiles above similarity threshold required for module search (e.g. 2); 
%                               if < 1 it is the fraction of profiles above threshold. 
%       N_RAND = integer number, e.g. 1000, number of randomization for each set
%       N_MAX = integer number, e.g. 4, number of bins on the number-of-profiles distribution for the randomizations
%
% OUTPUT: .mat and .txt files in the same folder as input files
%
% example call IMPACT_modules([1 1 1 1], [0.7  2 5000 4])
%
% by Giovanni Marsico and Angela Simeone 04/03/2014
% redistributed under The MathWorks, Inc. Software License Agreement

%% SETTING BASIC PARAMETER FOR ENVIRONMENT AND SEARCH

% environment parameters
do_seed_nodes = do_params(1);
do_module_search = do_params(2);
do_randomizations = do_params(3);
do_statistics = do_params(4);

% algorithm parameters
SIMIL_THR = alg_params(1);
MIN_NUMB_PROF = alg_params(2);
N_RAND = alg_params(3);
N_MAX = alg_params(4);
    
% set path and file names
[networkFileIN my_path] = uigetfile('*.txt', 'Load network interaction data','.');
networkFileIN = [my_path networkFileIN];
[ExperimentalDataFileIN my_path] = uigetfile('*.txt', 'Load phenotypic data',my_path);
ExperimentalDataFileIN  = [my_path ExperimentalDataFileIN];
[temp1 base_name] = fileparts(ExperimentalDataFileIN);

networkFile                 = [strtok(networkFileIN,'.') '.mat'];
ExperimentalDataFile        = [strtok(ExperimentalDataFileIN,'.') '.mat'];
SeedNodeFile                = [my_path base_name '_seed-file_prof-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.mat'];
VectorSeedNodeStartingFile  = [my_path base_name '_seed-list_prof-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.txt'];
ModuleResultsFile           = [my_path base_name '_res-search-file_prof-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.mat'];
OutputFileDistr             = [my_path base_name '_p-distr-file_prof-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.mat'];
OutputFileStats             = [my_path base_name '_res-stat-file_prof-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.mat'];
ExportTextFile              = [my_path base_name '_output-export-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100) '.txt'];
ListExportName              = [my_path base_name '_export_list-' num2str(MIN_NUMB_PROF) '_thr-' num2str(SIMIL_THR*100)];

if(do_seed_nodes || do_module_search)
    do_it = 1;
    if ( exist(networkFile, 'file') == 2)   % .mat file with the same name already exist
        button = questdlg('A network file for this input file has already been created. Do you want to overwrtie it?','Import network') ;
        if (strcmp(button,'No') == 1)
            do_it = 0;
        end
    end
    if (do_it)
        networkFile = convert_text_mat(1, 1, networkFileIN);
    end
    
    do_it = 1;
    if ( exist(ExperimentalDataFile, 'file') == 2)   % .mat file with the same name already exist
        button = questdlg('A phenotypic file for this input file has already been created. Do you want to overwrtie it?','Import profiles'); 
        if (strcmp(button,'No') == 1)
            do_it = 0;
        end
    end
    if (do_it)
        ExperimentalDataFile = convert_text_mat(1, 0, ExperimentalDataFileIN);
    end
end
%% FILE IMPORT, SEED NODE AND SEED PATTERN SELECTION

if (do_seed_nodes)
    [PotentialSeedNodeDef, SeedNodePatternsOriginal, SeedNodeProfilesOriginal] = RunTOPmoduleSeedSelection(networkFile, ExperimentalDataFile, 1, 0.8, SeedNodeFile, 2);
    vector=ones(size(PotentialSeedNodeDef));
    dlmwrite(VectorSeedNodeStartingFile,vector), clear vector;
end

%%  MODULE SEARCH
% input files required (string variables inside script)
%       netowrk file (e.g., 'HPRDinVivo_Intact_Kegg15092009.mat')
%       phenotypic data file (e.g., 'AllOligo40ProfilesEnsembl_oligSymb_Nozscore.mat')
%       seed node file (e.g., 'VectorFileCombinedNetSetseed08_numprof2_TESTJuly2012.txt')
%       seed pattern profile file (e.g., 'SeedNodeHPRDinVivo_numProf2_setNode08_TESTJuly2012.mat')
% output files 
%       module search result file (e.g., 'ResultsCombinedNetSetNodes2Prof_07simil_OCT2012.mat')
% parameters defined
%       SIMIL_THR = similarity theshold (e.g. 0.7);
%       MIN_NUMB_PROF = minimal number of required simialr profiles (e.g. 2);
if(do_module_search)
    module_search (networkFile, ExperimentalDataFile, VectorSeedNodeStartingFile, SeedNodeFile, ModuleResultsFile, SIMIL_THR, MIN_NUMB_PROF)
end

%% DISTRIBUTION CALCULATION 
% input files required (string variables inside script)
%       module search result file (e.g., 'ResultsCombinedNetSetNodes2Prof_08simil_OCT2012.mat')        
%       seed pattern profile file (e.g., 'SeedNodeHPRDinVivo_numProf2_setNode08_TESTJuly2012.mat')
% output files
%       probability distribution file for used seeds (e.g., 'P_seedPattern_ModeProfileOCT2012_1-100.mat')

if (do_randomizations)
    Probabilities_distribution_calculation(ModuleResultsFile,SeedNodeFile,OutputFileDistr, MIN_NUMB_PROF, SIMIL_THR, N_RAND, N_MAX);
end
%% MODULE STATISTICAL ASSESSMENT (module fusion, p-values, ROC curves)
% input files required (string variables inside script)
%       module search result file (e.g., 'ResultsCombinedNetSetNodes2Prof_08simil_OCT2012.mat')        
%       probability distribution file for used seeds (e.g., 'P_seedPattern_ModeProfileOCT2012_1-100.mat')
%       reference set (e.g., ReferenceSets_GOExperimental_GOAll_POP.mat)
% output files
%       p-value file (e.g., 'RESULTS_2prof_07simil_OCT2012.mat')

if (do_statistics)
    warning('off', 'MATLAB:nchoosek:LargeCoefficient');
    ModuleFusion_PvalueCalculation(ModuleResultsFile,OutputFileDistr,OutputFileStats,SIMIL_THR,N_MAX);
    warning('on', 'MATLAB:nchoosek:LargeCoefficient');
    th = 0.1;
    [list_all_network, list_all_modules, list_sig_modules] = PrintResultsToFile_Modules (OutputFileStats,ExportTextFile, th);
    name1 = [ListExportName '_allnet.txt'];
    name2 = [ListExportName '_allmod.txt'];
    name3 = [ListExportName '_allsignificant.txt'];
    cell2csv(name1,list_all_network);
    cell2csv(name2,list_all_modules);
    cell2csv(name3,list_sig_modules);
end