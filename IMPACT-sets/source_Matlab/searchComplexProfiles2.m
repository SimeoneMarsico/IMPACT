function [selectedProfiles,selectedProfilesNotSwap,referenceProfile,matrixProfilesAll,BestSize,NumGenesSelected,TotNumProf,selectGeneSymb,SelectedOligos]=searchComplexProfiles2(newComplexes,genesSymbols,oligoSymbol,THR_corr_user,simil_type)

%inizialization of all variable for storing  results
selectedProfiles=cell(1,length(newComplexes));
selectedProfilesNotSwap=cell(1,length(newComplexes));
referenceProfile=cell(1,length(newComplexes));

BestSize=zeros(1,length(newComplexes));
NumGenesSelected=zeros(1,length(newComplexes));
TotNumProf=zeros(1,length(newComplexes));
selectGeneSymb=cell(1,length(newComplexes));
SelectedOligos=cell(1,length(newComplexes));

%inizialization of calculation variables
matrixProfiles=[];
matrixProfilesAll={};
matrixOligos={};
matrixGenes=[];

h = waitbar(0,'searching enriched pattern...');
for i=1:length(newComplexes)
    if(mod(1000*i/length(newComplexes),10) == 0)
        waitbar(i/length(newComplexes));
    end
    if isempty(newComplexes{i})
        continue
    else
    complexSize=length(newComplexes{i});
    complexName=genesSymbols{i};
    complexOligos=oligoSymbol{i};
    
    %actual case
    for j=1:complexSize
        NprofGene=size(newComplexes{i}{j},1);
        matrixProfiles=[matrixProfiles;newComplexes{i}{j}];
        matrixGenes=[matrixGenes;ones(NprofGene,1)*j];
        matrixOligos=[matrixOligos;complexOligos{j}];
    end
 
    matrixProfiles2=matrixProfiles;
    matrixProfilesAll{i}=matrixProfiles;
    clear NprofGene j
    
    if(simil_type==0)
        CorrProf=corrcoef(matrixProfiles');
    else
        CorrProf=corrcoef_dist(matrixProfiles',0);
    end
    CorrProf=CorrProf-eye(size(CorrProf));
    
    for j=1:length(CorrProf)
        %indAboveThr=find(abs(CorrProf(j,:)) >= THR_corr_user);
        indAboveThr=find((CorrProf(j,:)) >= THR_corr_user);
        t(j)=length(indAboveThr);
		g(j)=length(unique([matrixGenes(indAboveThr);matrixGenes(j)]));
    end
    TotNumProf(i)=length(CorrProf);
    %searching of the row/col with the maximum number of similar profiles
    %selection of the row/col with the maximum number of similar profiles
    %showed by the maxum number of different genes
    % t: number of similar profiles (oligo based)
    % g: number of similar genes (gene based)
    flag = 1;
    while(flag) 
        [dir_]  = max(t); %max num of profiles
        inddir = find(t == dir_);
        flag = 0;
        for i1 = 1:length(inddir)
            all_corrs_in_compl = CorrProf(inddir(i1),:); 
           
            count = 0;
            pos = length(find(all_corrs_in_compl >= THR_corr_user));
            neg = length(find(all_corrs_in_compl <= -THR_corr_user));
            if(pos >= neg)
                inds = find(all_corrs_in_compl >= THR_corr_user);
            else
                inds = find(all_corrs_in_compl <= -THR_corr_user);
            end
            count = count + max(pos,neg);

            flag = flag || (count < dir_);
            if (count < dir_)
                t(inddir(i1)) = count;
                g(inddir(i1)) = length(unique([matrixGenes(inds);matrixGenes(inddir(i1))])); 
            end
        end
    end
    [dirG, inddirG] = max(g); %max num of genes
    
    
    %selection of the max number of profiles covering the max num of genes   
    inddir_=find(t==max(t));
    inddirG_=find(g==max(g));

    desideredRow=intersect(inddir_,inddirG_);
    if isempty(desideredRow)
        desideredRow=inddirG_;
    end
    if length(desideredRow)==1
        desideredRowFinal=desideredRow;
    end
    
    if length(desideredRow)>1 %case of multiple choises: select the one with the maximum of the sum(Corr)
        
        for pp=1:length(desideredRow)
            %colIndex=find(abs(CorrProf(desideredRow(pp),:))>THR_corr_user);
            %tempVect(pp)=sum(abs(CorrProf(desideredRow(pp),colIndex)));
            colIndex=find((CorrProf(desideredRow(pp),:))>THR_corr_user);
            tempVect(pp)=sum((CorrProf(desideredRow(pp),colIndex)));
            tempGenes(pp)=length(unique([matrixGenes(desideredRow(pp))' matrixGenes(colIndex)']));
        end
        
        [desideredRowMaxProf] = find(tempVect==max(tempVect));
        desideredRowMax=desideredRow(find(tempVect==max(tempVect))); % selection of the row with the maximum sum(corr)
        [desideredRowFinalGenes] = find(tempGenes==max(tempGenes));
        
        desideredRowtemp=desideredRow(intersect(desideredRowMaxProf,desideredRowFinalGenes));
        if length(desideredRowtemp)>1
            desideredRowFinal=setdiff(desideredRow,inddirG);
        elseif isempty(desideredRowtemp)
            desideredRowFinal=desideredRow(desideredRowFinalGenes(1));
        elseif length(desideredRowtemp)==1
            desideredRowFinal=desideredRowtemp;
        end
        clear colIndex tempVect tempGenes desideredRowtemp MaxProfdesideredRowMaxProf
    end
    
    %Code for:
    % - selection of profiles 
    % - swapping of anticorrelated profiles
    % - computation of the reference profile

        if dir_ >= 1 && dirG >=1

            indD=desideredRowFinal;
            
            if length(indD)>1
                indD=indD(1); % this has to be checked!!!there is a problem! some times indD contains multiple indeces that means that there are computations of the same couple of profiles multiple times!

            end
            %corr_set=find(abs(CorrProf(indD,:))>THR_corr_user); %indeces vector
            corr_set=find((CorrProf(indD,:))>THR_corr_user); %indeces vector
            indprofilesProva=[indD corr_set];

            % code for selecting the profiles that nees to be swapped
            corr_sign=sign(CorrProf(indD,corr_set));

            %the following vectors contain indeces
            pos_corr=find(corr_sign==1);
            neg_corr=find(corr_sign==-1);

            if length(pos_corr)>length(neg_corr)
                indToSwap=corr_set(neg_corr);
            elseif length(pos_corr)<length(neg_corr)
                indToSwap=[indD corr_set(pos_corr)];
            elseif length(pos_corr)==length(neg_corr)
                sumPos=sum(abs(CorrProf(indD,corr_set(pos_corr))));
                sumNeg=sum(abs(CorrProf(indD,corr_set(neg_corr))));
                if sumPos>sumNeg
                    indToSwap=corr_set(neg_corr);
                elseif sumPos<sumNeg
                    indToSwap=[indD corr_set(pos_corr)];
                end
            end

            matrixProfiles2(indToSwap,:)=matrixProfiles2(indToSwap,:)*-1;

            clear maxt indt 

            selectedProfiles{i}=matrixProfiles2(indprofilesProva,:);
            selectedProfilesNotSwap{i}=matrixProfiles(indprofilesProva,:);
            referenceProfile{i}=median(matrixProfiles2(indprofilesProva,:));
            if(simil_type==0)
                tempBestSize=corrcoef([referenceProfile{i};selectedProfiles{i}]');
            else
                tempBestSize=corrcoef_dist([referenceProfile{i};selectedProfiles{i}]',0);
            end
            BestSize(i)=length(indprofilesProva);

            SelectedOligos{i}=matrixOligos(indprofilesProva);
            clear tempBestSize

            indexGenes=unique(matrixGenes(indprofilesProva));
            if (~iscell(complexName))
                complexName = cellstr(complexName);
            end
            selectGeneSymb{i}=complexName(indexGenes);
            ll(i)=length(indexGenes);
            percentage=ll/length(newComplexes{i})*100;

        clear indToSwap pos_corr neg_corr  sumPos sumNeg ambig_corr corr_set corr_sign percentage indexGenes

        end
    
    end % close if isempty loop
    
    clear dir inddir
    clear tempVect pp desideredRow desideredRowFinal
    clear NprofGene matrixGenes  CorrProf matrixProfiles k t maxdir indD indprofilesProva percentage g 
    matrixProfiles=[];
    matrixGenes=[];
    matrixOligos={};
    
end

NumGenesSelected=ll; clear ll
clear i j matrixProfiles matrixGenes matrixProfiles2
if(exist('h'))
    close(h);
end


        