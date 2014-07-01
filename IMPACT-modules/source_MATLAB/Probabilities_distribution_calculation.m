function Probabilities_distribution_calculation(ModuleResultsFile,SeedNodeFile,OutputFile, n_min, simil_thr, n_rand, N_max, simil_type)
    
    load(ModuleResultsFile);
    load(SeedNodeFile);
    h = waitbar(0,'calculating distribution probability');
    P_SeedPatternUsed = cell(1,length(SeedPatternUsed));
    for i=1:length(SeedPatternUsed)
        waitbar(i/length(SeedPatternUsed));
        P_SeedPatternUsed{i} = threshold_frequency(AllPROFnet,SeedPatternUsed(i,:), n_rand, simil_thr, n_min, N_max, simil_type);
    end
    if (exist('h'))
        close(h);
    end
    save(OutputFile,'P_SeedPatternUsed','SeedPatternUsed');
end