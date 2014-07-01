function [flag_input] = import_map_data(interactions_data_file,multiparam_data_file, imported_mapped_dataFileName, multiparametricData_matFileName)
% IMPORT DATA FROM TEXT FILE
% function [flag_input,imported_data] = import_map_data(interaction_data_file,multiparam_data_file)
% This function requires two text files as inputs.
% the function reads the phenotypic data and procedes with the mapping on
% the interaction data.
% The input files must be text files containing info like

% newComplexes,newComplexesID,newComplexesSymbol,SIMIL_THR

fid = fopen(multiparam_data_file,'rt');
a_line = fgetl(fid);
fclose(fid)
[a b] = strtok(a_line);
count = 1;
while(~isempty(b))
    [a b] = strtok(strtrim(b));
    count = count + 1;
end
numParamToUse = count - 2;

fid_interactions = fopen(interactions_data_file,'rt','n','UTF-8');
if fid_interactions < 0
    disp('**ERROR ** --- Invalid interactions data filename')
    return
end
disp (' ')
disp (' ')
disp ('Import set information -- started ')
[Set_info] = textscan(fid_interactions, '%s %s %s %s %s', 'delimiter', '\t');
DatabaseID = Set_info{1}; % it is assumed there is a header for each column
Set_Name = Set_info{2};
Set_Components = Set_info{3};
Set_Organism = Set_info{4};
Set_dbid = Set_info{5};

fclose(fid_interactions); clear fid_interactions Set_info
h = waitbar(0,'importing sets...');
%Set_Components_Details = cell(1,length(Set_Components));
for i = 2:length(Set_Components)
    if(mod(1000*i/length(Set_Components),10) == 0)
        waitbar(i/length(Set_Components));
    end
    ind_Commas = strfind(Set_Components{i},',');
    indices = [0 ind_Commas length(Set_Components{i})];
    
    if isempty(ind_Commas)
        Set_Components_Details{i} = Set_Components{i};
    else
        tempComponent = Set_Components{i};
        for j = 1:length(indices)-1
            [Set_Components_Details{i}{j} tempComponent]= strtok(tempComponent,',');
        end
        clear ind_Commas indices tempComponent
    end
end
if(exist('h'))
    close(h);
end
        

disp ('... terminated succesfully')
disp (' ')
disp (' ')
disp ('Import multiparametric/phenotypic information -- started ')
fid_phentypic = fopen(multiparam_data_file,'r');
if fid_phentypic < 0
    disp('**ERROR ** --- Invalid phenotypic data filename')
    return
end


Params_to_read = repmat(' %f', 1, numParamToUse);
Format_to_read = ['%s %s' Params_to_read];
Data_prof_symb = textscan(fid_phentypic,Format_to_read, 'delimiter', '\t');
DataID = unique(Data_prof_symb{1});
Data_Prof_temp =  Data_prof_symb(3:end);
Data_Prof = cell2mat(Data_Prof_temp); % matrix of profiles
fclose(fid_phentypic);
clear Data_Prof_temp Params_to_read Format_to_read fid_phentypic
disp ('... terminated succesfully.')
disp (' ')
disp ('Mapping step (mapping multiparametric information onto set data -- started')
DataIDProfiles = cell(size(DataID));
DataIDsymbols =  cell(size(DataID));
sizeID = length(DataID{1});
m=1;
h = waitbar(0,'mapping pehnotypic data onto sets...');
for i=1:length(DataID)
    if(mod(1000*i/length(DataID),10) == 0)
        waitbar(i/length(DataID));
    end
    currentID = DataID{i};
    currentIDoccurances = strncmp(currentID, Data_prof_symb{1},sizeID);
    indCurrentID = find(currentIDoccurances);
    currentID_profiles= Data_Prof (indCurrentID,:);
    DataIDProfiles{i} = currentID_profiles;
    DataIDsymbols{i} = Data_prof_symb{2}(indCurrentID);
    clear currentID indCurrentID currentID_profiles currentIDoccurances
    m=m+1;
end
close(h);
clear m

allOligoProfID = DataID;
allOLIGOName = DataIDsymbols;
allOligoProf = DataIDProfiles;



% map phenotypic data onto set information (interaction data)
% newComplexes,newComplexesID,newComplexesSymbol
clear m
m=1;

h = waitbar(0,'creating output files...');

for i=1:length(Set_Name)
    if(mod(1000*i/length(Set_Name),10) == 0)
        waitbar(i/length(Set_Name));
    end
    if not(iscell(Set_Components_Details{i}))
        currentComponent = Set_Components_Details{i};
        indCurrentComponent = find(strncmp(currentComponent,DataID,length(currentComponent)));
        if not(isempty(indCurrentComponent))
            SetProfile{m} = DataIDProfiles (indCurrentComponent);
            SetSymbol{m} = DataIDsymbols (indCurrentComponent);
            SetName{m} = Set_Name{i};
            SetID{m} = Set_Components_Details{i};
            Setdbid{m} = Set_dbid{i};
            SetOrganisms{m} = Set_Organism{i};
            m=m+1;

        end
    else 
        SetName{m} = Set_Name{i};
        SetID{m} = Set_Components_Details{i};
        count_j = 1;
        for j=1:length(Set_Components_Details{i});
            currentComponent = Set_Components_Details{i}{j};
            indCurrentComponent = find(strncmp(currentComponent,DataID,length(currentComponent)));
            if(~isempty(indCurrentComponent))
                SetProfile{m}{count_j} = DataIDProfiles{indCurrentComponent};
                SetSymbol{m}{count_j} = DataIDsymbols{indCurrentComponent};
                Setdbid{m} = Set_dbid{i};
                SetOrganisms{m} = Set_Organism{i};
                count_j = count_j + 1;
            end
        end
        m=m+1;
    end
end
close(h);
newComplexes = SetProfile;
newComplexesID = SetID;
newComplexesName = SetName;
newComplexesSymbol = SetSymbol;
newComplexesIDdb = Setdbid;
newComplexesOrganism = SetOrganisms;
TotalSetID=[];
for i=1:length(newComplexesID)
    a = newComplexesID{i};
    if not(iscell(a))
        a = {a};
    end
    TotalSetID = [TotalSetID a];
end
    

allProteinSetUnique= unique(TotalSetID); clear TotalSedID
[allProteinSetUniqueID indUnique] = intersect(allProteinSetUnique,allOligoProfID);

allProteinSetSymbolsUnique = allOLIGOName(indUnique); 
allProteinSetUnique = allOligoProf(indUnique);


disp (' ')
disp ('.. terminated succesfully.')

% imported_data.profile = SetProfile;
% imported_data.symbols = SetSymbol;
% imported_data.Name = SetName;
% imported_data.ID = SetID;
flag_input = 1;
disp (' ')
disp (' ')
disp (['Variables stored in output file:' imported_mapped_dataFileName])
disp (' ')
save(imported_mapped_dataFileName,'newComplexes*','allProteinSet*','allOligoProf');
save(multiparametricData_matFileName,'allOligoProfID','allOLIGOName','allOligoProf');

    

