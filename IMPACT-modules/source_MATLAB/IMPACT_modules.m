function IMPACT_modules
% loaded input file parameters should look like:
% do_seed_nodes,do_module_search,do_randomizations,do_statistics,SIMIL_THR,MIN_NUMB_PROF,N_RAND,N_MAX
% e.g. 1,1,1,1,0.7,2,5000,4
% where do_seed_nodes     = binary, [0,1], do the import (convert input text files into matlab environment files)
%       do_module_search  = binary, [0,1], do the search procedure on sets
%       do_randomizations = binary, [0,1], do the randomizations
%       do_statistics     = binary, [0,1], calculate p-values and export results based on the randomization results
%       simil_type        = binary, [0,1], type of similarity measure: if 0 is Pearson correlation coefficient, if 1 is inverse of Euclidean distance
%       SIMIL_THR = number, [0-1], similarity threshold (Pearson correlation coefficient) for set search and randomizations
%       MIN_NUMB_PROF = number, if >= 1 it is the number of profiles above similarity threshold required for module search (e.g. 2); 
%                               if < 1 it is the fraction of profiles above threshold. 
%       N_RAND = integer number, e.g. 1000, number of randomization for each set
%       N_MAX = integer number, e.g. 4, number of bins on the number-of-profiles distribution for the randomizations
%
% OUTPUT: .mat and .txt files in the same folder as input files

%% SETTING BASIC PARAMETER FOR ENVIRONMENT AND SEARCH

[param_file param_file_path] = uigetfile('*.txt', 'Load paramters file','.');
a = dlmread([param_file_path param_file]);
do_seed_nodes = a(1)
do_module_search = a(2)
do_randomizations = a(3)
do_statistics = a(4)
simil_type = a(5)
SIMIL_THR = a(6)
MIN_NUMB_PROF = a(7)
N_RAND = a(8)
N_MAX = a(9)
 
% set path and file names
[networkFileIN my_path] = uigetfile('*.txt', 'Load network interaction data','.');
networkFileIN = [my_path networkFileIN];
[ExperimentalDataFileIN my_path] = uigetfile('*.txt', 'Load phenotypic data',my_path);
ExperimentalDataFileIN  = [my_path ExperimentalDataFileIN];
[temp1 base_name] = fileparts(ExperimentalDataFileIN);

%networkFile                 = [strtok(networkFileIN,'.') '.mat'];
%ExperimentalDataFile        = [strtok(ExperimentalDataFileIN,'.') '.mat'];
all_dots = strfind(networkFileIN,'.');
networkFile = [networkFileIN(1:all_dots(end)-1) '.mat'];
all_dots = strfind(ExperimentalDataFileIN,'.');
ExperimentalDataFile = [ExperimentalDataFileIN(1:all_dots(end)-1) '.mat'];

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

T_s = 3;
k_s = 1;
if (do_seed_nodes)
    %[PotentialSeedNodeDef, SeedNodePatternsOriginal, SeedNodeProfilesOriginal] = RunTOPmoduleSeedSelection(networkFile, ExperimentalDataFile, 1, 0.8, SeedNodeFile, 3);
    [PotentialSeedNodeDef, SeedNodePatternsOriginal, SeedNodeProfilesOriginal] = RunTOPmoduleSeedSelection(networkFile, ExperimentalDataFile, 1, T_s, SeedNodeFile, k_s, simil_type); % crispr was 3 and 3; autophagy was 4...1; depletion was 5 and 1
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
    module_search (networkFile, ExperimentalDataFile, VectorSeedNodeStartingFile, SeedNodeFile, ModuleResultsFile, SIMIL_THR, MIN_NUMB_PROF,simil_type)
end

%% DISTRIBUTION CALCULATION 
% input files required (string variables inside script)
%       module search result file (e.g., 'ResultsCombinedNetSetNodes2Prof_08simil_OCT2012.mat')        
%       seed pattern profile file (e.g., 'SeedNodeHPRDinVivo_numProf2_setNode08_TESTJuly2012.mat')
% output files
%       probability distribution file for used seeds (e.g., 'P_seedPattern_ModeProfileOCT2012_1-100.mat')

if (do_randomizations)
    Probabilities_distribution_calculation(ModuleResultsFile,SeedNodeFile,OutputFileDistr, MIN_NUMB_PROF, SIMIL_THR, N_RAND, N_MAX, simil_type);
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
    ModuleFusion_PvalueCalculation(ModuleResultsFile,OutputFileDistr,OutputFileStats,SIMIL_THR,N_MAX, simil_type);
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