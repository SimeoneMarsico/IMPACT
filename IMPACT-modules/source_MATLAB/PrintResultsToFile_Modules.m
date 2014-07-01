function [list_all_network, list_all_modules, list_sig_modules] = PrintResultsToFile_Modules (MAT_ModuleResultsFileMat,TXT_outputfile, th)

load(MAT_ModuleResultsFileMat)
fid = fopen(TXT_outputfile,'w');

list_all_network = AllGenesNet_Screen;
list_all_modules = {};
list_sig_modules = {};

fprintf(fid, 'NUM_MODULE\tp-value\t#tot_genes\t#tot_profiles\tGene_ID\tProfiles_Symbols\n');
for i=1:length(moduleSet)
    if ~isempty(moduleSet{i})
        pval = modulePvalue(i);
        tot_genes = length(moduleSet{i});
        tot_profiles = size(moduleProfile{i},1);
        GeneID = sprintf('%s,',AllGenesNet_Screen{moduleSet{i}});
        list_all_modules = [list_all_modules;AllGenesNet_Screen(moduleSet{i})];
        if (pval <= th)
            list_sig_modules = [list_sig_modules; AllGenesNet_Screen(moduleSet{i})];
        end
        Profiles_Symbols = sprintf('%s,',moduleOligoSymb{i}{:});
        fprintf(fid,'%i\t%f\t%i\t%i\t%s\t%s\n',i,pval,tot_genes,tot_profiles,GeneID,Profiles_Symbols);
        clear pval tot_genes tot_profiles GeneID Profiles_Symbols
    end
end

fclose(fid)
clear fid
list_all_modules = unique(list_all_modules);
list_sig_modules = unique(list_sig_modules);
        
        