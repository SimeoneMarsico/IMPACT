function IMPACT_sets(do_params, alg_params)
% INPUT:
% do_params = [do_import_mapping do_set_search do_randomizations do_statistics]
% where do_import_mapping = binary, [0,1], do the seed node finding
%       do_set_search = binary, [0,1], do the search procedure on network
%       do_randomizations = binary, [0,1], do the randomizations
%       do_statistics = binary, [0,1], calculate p-values and export results based on the randomization results
% alg_params = [SIMIL_THR MIN_NUMB_PROF N_RAND N_MAX]
% where SIMIL_THR = number, [0-1], similarity threshold (Pearson correlation coefficient) for set search and randomizations
%       MIN_NUMB_PROF = not used, kept only for compatibility with IMPACT-modules.m (IMPACT-modules); it must be passed to the function (it may be used in future implementations);
%       N_RAND = integer number, e.g. 1000, number of randomization for each set
%       N_MAX = integer number, e.g. 4, number of bins on the number-of-profiles distribution for the randomizations
%
% OUTPUT: .mat and .txt files in the same folder as input files
%
% example call IMPACT_sets([1 1 1 1], [0.75 1 5000 4])
%
% by Giovanni Marsico and Angela Simeone 04/03/2014
% redistributed under The MathWorks, Inc. Software License Agreement


%% SETTING BASIC PARAMETER FOR ENVIRONMENT AND SEARCH

% environment parameters
do_import_mapping = do_params(1);
do_set_search = do_params(2);
do_randomizations = do_params(3);
do_statistics = do_params(4);

% algorithm parameters
SIMIL_THR = alg_params(1);
MIN_NUMB_PROF = alg_params(2);
N_RAND = alg_params(3);
N_MAX = alg_params(4);
    
[set_interaction_file my_path] = uigetfile('*.txt', 'Load set interaction data','.');
set_interaction_file = [my_path set_interaction_file];
[experimental_pheotypic_file my_path] = uigetfile('*.txt', 'Load phenotypic data',my_path);
experimental_pheotypic_file  = [my_path experimental_pheotypic_file];
my_slash = my_path(end);
[temp1 base_name] = fileparts(experimental_pheotypic_file);

output_file_pheno_mapped = [my_path base_name '_pheno-mapped.mat'];
output_file_pheno = [my_path 'MultiparametricDataAllprofile-08_ImportedTEST.mat'];  % useless remove later on
RandomizationsFolder = [my_path 'RandomizationsFolder' my_slash];
if not(isdir(RandomizationsFolder))
    mkdir(RandomizationsFolder)
end
OutputRandFile   = [RandomizationsFolder 'Randomizations_thr-' num2str(100*SIMIL_THR) '_n-rand-' num2str(N_RAND) '.mat'];
OutputSearchFile = [my_path base_name 'res-search_thr-' num2str(100*SIMIL_THR) '.mat'];
OutputResultFile = [my_path base_name '_res-final_thr-' num2str(100*SIMIL_THR) '.mat'];
ExportTextFile   = [my_path base_name '_output_export-' num2str(100*SIMIL_THR) '.txt'];
ListExportName   = [my_path base_name '_export_list-'   num2str(SIMIL_THR*100)];

%% import data from txt file and create mat file to use for next step analysis
if(do_import_mapping)
    do_it = 1;
    if ( exist(output_file_pheno_mapped, 'file') == 2)   % .mat file with the same name already exist
        button = questdlg('A mapping file for this input file has already been created. Do you want to overwrtie it?','Import/Map profiles') ;
        if (strcmp(button,'No') == 1)
            do_it = 0;
        end
    end
    if (do_it)
        [flag_input] = import_map_data(set_interaction_file,experimental_pheotypic_file, output_file_pheno_mapped, output_file_pheno);
    end
end
% Load phenotypic data mapped on the interaction data (sets, e.g. protein complexes from CORUM)
if (do_set_search || do_randomizations || do_statistics)
    load(output_file_pheno_mapped)
end
%% Set-base analysis: search enriched pattern in complexes 
if(do_set_search)
    [Realselected,RealselectedPnotSwap,RealreferenceP,RealmatrixP,RealMaxNUmProf,RealNumGenesSelected,RealTotNumProf,RealSelGenesID,RealSelectedOligos]=searchComplexProfiles2(newComplexes,newComplexesID,newComplexesSymbol,SIMIL_THR);
    save (OutputSearchFile, 'Real*');
else
    load(OutputSearchFile);
end

%% RANDOMIZATION PROCEDURE
% this steps allow to estimate the statistical significance by computing a
% pvalue using a permutation strategy
% 1 - run searching procedure over a randomly assembled set
if (do_randomizations)
    
    % all screens profiles
    numProfilesPerGene_unique = zeros(1,length(allOligoProf));
    for i=1:length(allOligoProf)
        numProfilesPerGene_unique(i) = size(allOligoProf{i},1);
    end
    % creation of 
    % Bins : set of bins containing gene profiles according to their number of profiles 
    % vectRandTemp : total count of genes in each bin.
    
    temp = sort(numProfilesPerGene_unique);
    for i = 1: N_MAX
        N(i) = temp(i * floor(length(temp)   /(N_MAX+1)));
    end
    N = unique(N);
    N_MAX = min(N_MAX, length(N));
    for i = 1 : N_MAX
        if i == 1
            Bins{i}=allOligoProf(numProfilesPerGene_unique<=N(i));
            vectRandTemp(i,:)=[length(find(numProfilesPerGene_unique<=N(i))) 1];
        elseif i < N_MAX
            Bins{i}=allOligoProf(numProfilesPerGene_unique > N(i-1) & numProfilesPerGene_unique <= N(i));
            vectRandTemp(i,:)=[length(find(numProfilesPerGene_unique > N(i-1) & numProfilesPerGene_unique <= N(i))) i];
        else
            Bins{i}=allOligoProf(numProfilesPerGene_unique>=N(i));
            vectRandTemp(N_MAX,:)=[length(find(numProfilesPerGene_unique>=N(i))) N_MAX];
        end
    end
    
    SetProfiles = newComplexes;
    BinnedSetinfo_PhenotypicDataFile = [my_path 'BinnedDatasetForRandomizations.mat'];
    save(BinnedSetinfo_PhenotypicDataFile, 'SetProfiles', 'Bins', 'vectRandTemp');

    RunRandomizationWithBinningSet(BinnedSetinfo_PhenotypicDataFile,OutputRandFile,N_RAND,SIMIL_THR, N);

end

%% STATISTICS CALCULATION BASED ON RANDOMIZATIONS
if (do_statistics)

    % Load results of randomizations steps
    Randomizations_selectedGenes=[];
    Randomizations_selectedProfiles=[];
    tempFilenames=dir(RandomizationsFolder);
    tempFilenames2={tempFilenames.name};
    FilenamesInd=strmatch('Randomiz',tempFilenames2);
    Filenames=tempFilenames2(FilenamesInd);
    clear tempFilenames tempFilenames2 FilenamesInd

    for i=1:length(Filenames)
        load([RandomizationsFolder Filenames{i}])
        Randomizations_selectedGenes=[Randomizations_selectedGenes;Rand_SelectedGenes];
        Randomizations_selectedProfiles=[Randomizations_selectedProfiles;Rand_selectedProfiles];
    end

    % P-value computations
    pvalueGenes = pvalue_fromFrequenciesVersionNovember09(RealreferenceP,RealNumGenesSelected,Randomizations_selectedGenes);
    pvalueProfiles = pvalue_fromFrequenciesVersionNovember09(RealreferenceP,RealMaxNUmProf,Randomizations_selectedProfiles);

    % preparation and saving output files
    selectedGenesAnalysis = {};
    selectedPvalues = [];
    for i=1:length(RealSelGenesID)
        if ~isempty(RealreferenceP{i})
            if (~iscell(RealSelGenesID{i}))
                if (length(RealSelGenesID{i} == 1))
                     RealSelGenesID{i} = cellstr(newComplexesID{i});
                else
                    RealSelGenesID{i} = cellstr(RealSelGenesID{i});
                end
            end
            for j=1:length(RealSelGenesID{i})
                if ~isempty(RealSelGenesID{i}{j})
                    selectedGenesAnalysis = [selectedGenesAnalysis; RealSelGenesID{i}{j}];
                    selectedPvalues = [selectedPvalues; pvalueGenes(i)'];
                end
            end
        end
        clear ngenes values
    end
    
    save(OutputResultFile, 'pvalue*', 'Real*', 'newComplexes*', 'allProteinSet*');
    th = 0.1;
    [list_all_sets, list_sel_sets, list_sig_sets] = PrintResultsToFile_Complexes (OutputResultFile,ExportTextFile, th);
    name1 = [ListExportName '_all-sets.txt'];
    name2 = [ListExportName '_sel-sets.txt'];
    name3 = [ListExportName '_sign-sets.txt'];
    cell2csv(name1,list_all_sets);
    cell2csv(name2,list_sel_sets);
    cell2csv(name3,list_sig_sets);
end