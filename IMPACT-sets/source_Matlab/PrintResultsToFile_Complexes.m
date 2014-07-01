function [list_all_sets, list_sel_sets, list_sig_sets] = PrintResultsToFile_Complexes (MAT_ComplexResultsFileMat,TXT_outputfile, th)

load(MAT_ComplexResultsFileMat)

list_all_sets = allProteinSetUniqueID;
list_sel_sets = {};
list_sig_sets = {};

fid = fopen(TXT_outputfile,'w');
fprintf(fid, 'NUM_SET\tSet_Name\tp-value(Gene)\tp-value(Profiles)\t#tot_genes_selected\t#tot_profiles_selected\tGene_ID(Selected)\tProfiles_Symbols(Selected)\n');
for i=1:length(RealMaxNUmProf)
    if RealMaxNUmProf(i)>0
        setName = newComplexesName{i};
        pvalueG = pvalueGenes(i);
        pvalueP = pvalueProfiles(i);
        tot_genes_sel = RealNumGenesSelected(i);
        tot_profiles_sel = RealMaxNUmProf(i);
        GeneID = sprintf('%s,',RealSelGenesID{i}{:});
        list_sel_sets = [list_sel_sets, RealSelGenesID{i}];
        if (pvalueG <= th)
            list_sig_sets = [list_sig_sets, RealSelGenesID{i}];
        end
        Profiles_Symbols = sprintf('%s,',RealSelectedOligos{i}{:});
        fprintf(fid,'%i\t%s\t%f\t%f\t%i\t%i\t%s\t%s\n',i,newComplexesName{i}, pvalueG, pvalueP, tot_genes_sel, tot_profiles_sel, GeneID, Profiles_Symbols);
        clear pval tot_genes tot_profiles GeneID Profiles_Symbols
    end
end

fclose(fid)
clear fid
list_sel_sets = unique(list_sel_sets);
list_sig_sets = unique(list_sig_sets);
if (size(list_all_sets,2) > size(list_all_sets,1))
    list_all_sets = list_all_sets';
end
if (size(list_sel_sets,2) > size(list_sel_sets,1))
    list_sel_sets = list_sel_sets';
end
if (size(list_sig_sets,2) > size(list_sig_sets,1))
    list_sig_sets = list_sig_sets';
end