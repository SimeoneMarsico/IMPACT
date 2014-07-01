function out = convert_text_mat(direction, type, filename)
% direction: 1 is text to mat, 0 is mat to text
% type: 1 is network type file, 0 is profiles
% filename: string of full file name including path
% out = output file conatining the full name of the converted file
% example call: out = convert_text_mat(1, 1, 1, 'my_path\my_filename');

    if (direction)
        if(type)
            h = waitbar(0,'converting network file from .txt to .mat');
            networkFile                 = filename;
            [HiConfNet(:,1), HiConfNet(:,2)]=textread(networkFile, '%s %s', 'delimiter', '\t');
            %alternatively, if textread will be removed from Matlab functions set
            % fid=fopen(networkFile,'rt','n','UTF-8')
            % HiConfNet = textscan(fid, '%s\t%s');
            % fclose(fid)
            % clear fid
            AllGenesNet = union(HiConfNet(:,1), HiConfNet(:,2));
            AdjMatrixNet = zeros(length(AllGenesNet));
            ht = java.util.Hashtable;          
            for i = 1:length(AllGenesNet)                                                     
                ht.put(AllGenesNet{i},i);          
            end
            AdjMatrixNet=zeros(size(AdjMatrixNet));
            for i=1:length(HiConfNet)
                if (mod(1000 * i/length(HiConfNet),10) == 0)
                    waitbar(i/length(HiConfNet));
                end
                temp1 = HiConfNet{i,1};
                temp2 = HiConfNet{i,2};
                ind1 = ht.get(temp1);
                ind2 = ht.get(temp2);
                AdjMatrixNet(ind1,ind2) = 1;
                AdjMatrixNet(ind2,ind1) = 1;
                clear temp1 temp2 ind1 ind2
            end
            AdjMatrixNet = sparse(AdjMatrixNet); % ok
            clear filename i ht
            [base_path base_name] = fileparts(networkFile);
            out = fullfile(base_path, [base_name '.mat']);
            save(out,'AdjMatrixNet','AllGenesNet');
            close(h); 
        else
            ExperimentalDataFile = filename;
            fid = fopen(ExperimentalDataFile,'rt');
            a_line = fgetl(fid);
            fclose(fid)
            [a b] = strtok(a_line);
            count = 1;
            while(~isempty(b))
                [a b] = strtok(strtrim(b));
                count = count + 1;
            end
            numParamToUse = count - 2;
            fid = fopen(ExperimentalDataFile,'rt');
            Params_to_read = repmat(' %f', 1, numParamToUse);clear numParamToUse
            Format_to_read = ['%s %s' Params_to_read]; clear Params_to_read
            Data_prof_symb = textscan(fid,Format_to_read, 'delimiter', '\t'); clear Format_to_read
            % multiple: 1 is multiple, 0 is single profile per gene
            multiple = 1;
            a = unique(Data_prof_symb{1});
            if (length(a) == length (Data_prof_symb{1}))
                multiple = 0;
            end
            if (~multiple)
                h = waitbar(0,'converting phenotypic file from .txt to .mat');
                % transform .txt to .mat - single profile per gene
                %path = 'C:\Users\marsico\Desktop\Angela_experiment_summary\review_comparisons\Code_NetworkAnalysis_1_07\';
                %ExperimentalDataFile        = 'AllOligo40ProfilesEnsembl_oligSymb_Nozscore.txt';
                
                fclose(fid)
                allOligoEns = Data_prof_symb{1};
                allOLIGOName = Data_prof_symb{2}';
                for i = 1:length(allOLIGOName)
                    allOLIGOName{i} = {allOLIGOName{i}};
                end
                Data_Prof = cell2mat(Data_prof_symb(3:end));
                allOligoProf = cell(1,length(allOligoEns));
                for i=1:length(allOligoEns)
                    allOligoProf{i} = Data_Prof(i,:);
                end
                clear fid Data_prof_symb Params_to_read Format_to_read Data_Prof
                close(h);
            else
                % transform .txt to .mat - multiple profiles per gene
                h = waitbar(0,'converting phenotypic file from .txt to .mat');
                DataID = unique(Data_prof_symb{1});
                Data_Prof_temp =  Data_prof_symb(3:end);
                Data_Prof = cell2mat(Data_Prof_temp); % matrix of profiles
                fclose(fid);

                sizeID = length(DataID{1});
                allOligoEns = DataID;  %#ok<NASGU>
                for i=1:length(DataID)
                    if(mod(i,100) == 0)
                        i
                    end
                    currentID = DataID{i};
                    currentIDoccurances = strncmp(currentID, Data_prof_symb{1},sizeID);
                    indCurrentID = find(currentIDoccurances);
                    currentID_profiles= Data_Prof (indCurrentID,:);
                    allOligoProf{i} = currentID_profiles;
                    allOLIGOName{i} = Data_prof_symb{2}(indCurrentID);
                    clear currentID indCurrentID currentID_profiles currentIDoccurances
                end
                clear sizeID
                
                close(h);
            end
            [base_path base_name] = fileparts(ExperimentalDataFile);
            out = fullfile(base_path, [base_name '.mat']);
            save(out,'allOligoEns','allOLIGOName', 'allOligoProf');
        end    
    else
        if(type)
            % transform .mat in text
            h = waitbar(0,'converting network file from .mat to .txt');
            networkFile                 = filename;
            load(networkFile);
            fid = fopen([strtok(networkFile,'/.') '.txt'], 'w','n','UTF-8');
            for i = 1:size(AdjMatrixNet,1)
                if (mod(1000*i/size(AdjMatrixNet,1),10) == 0)
                    waitbar(i/size(AdjMatrixNet,1));
                end
                for j = 1:size(AdjMatrixNet,2)
                    if (AdjMatrixNet(i,j))
                        fprintf(fid,'%s\t%s\n',AllGenesNet{i}, AllGenesNet{j});
                    end
                end
            end
            fclose(fid);
            clear fid;
            close(h);
            out = [strtok(networkFile,'/.') '.txt'];
        else
            ExperimentalDataFile = filename;
            load(ExperimentalDataFile);
            fid=fopen([strtok(ExperimentalDataFile,'/.') '.txt'],'w','n','UTF-8');
            h = waitbar(0,'converting phenotypic file from .mat to .txt');
            for i=1:length(allOligoEns)
                if (mod(1000*i/length(allOligoEns),10) == 0)
                    waitbar(i/length(allOligoEns));
                end
                for i1 = 1:length(allOLIGOName{i})
                    fprintf(fid, '%s\t%s\t', allOligoEns{i}, allOLIGOName{i}{i1});
                    temp = allOligoProf{i}(i1,:);
                    for j=1:length(temp)
                        fprintf(fid,'%s\t',num2str(temp(j)));
                    end
                    fprintf(fid,'\n');
                end
            end 
            fclose(fid)
            close(h);
            out = [strtok(ExperimentalDataFile,'/.') '.txt'];
        end
    end
end