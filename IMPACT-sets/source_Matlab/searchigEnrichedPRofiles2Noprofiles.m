function [bestSize,numSelGenes]=searchigEnrichedPRofiles2Noprofiles(complex,ThrCorr,simil_type)

matrixProfiles=[];
matrixGenes=[];

complexSize=length(complex);
currComplex=complex;
%actual case
for j=1:complexSize
    NprofGene=size(currComplex{j},1);
    %matrixProfiles=[matrixProfiles;currComplex{j}(:,1:20)];
    matrixProfiles=[matrixProfiles;currComplex{j}];
    matrixGenes=[matrixGenes;ones(NprofGene,1)*j];
end

matrixProfiles2=matrixProfiles;
clear NprofGene j
if(simil_type == 0)
    CorrProf=corrcoef(matrixProfiles');
else
    CorrProf=corrcoef_dist(matrixProfiles',0);
end
CorrProf=CorrProf-eye(size(CorrProf)); %substract the diagonal that contains 1

%preventiveCheck=find(abs(CorrProf) > ThrCorr);
preventiveCheck=find((CorrProf) > ThrCorr);
if preventiveCheck
    %if lengthpreventiveCheck
    for j=1:length(CorrProf)
        %indAboveThr=find(abs(CorrProf(j,:)) >= ThrCorr);
        indAboveThr=find((CorrProf(j,:)) >= ThrCorr);
        t(j)=length(indAboveThr); % number of profiles (along the rows) above threshold
		g(j)=length(unique(matrixGenes(indAboveThr))+1); % number of genes above threshold
    end
    %TotNumProf(i)=length(CorrProf);
    %searching of the row/col with the maximum number of similar profiles
    %selection of the row/col with the maximum number of smilar profiles
    %showed by the maxum number of different genes
    %[dir,inddir]=max(t);
    flag = 1;
    while(flag) 
        [dir_]  = max(t); %max num of profiles
        inddir = find(t == dir_);
        flag = 0;
        for i1 = 1:length(inddir)
            all_corrs_in_compl = CorrProf(inddir(i1),:);
           
            count = 0;
            pos = length(find(all_corrs_in_compl >= ThrCorr));
            neg = length(find(all_corrs_in_compl <= -ThrCorr));
            if(pos >= neg)
                inds = find(all_corrs_in_compl >= ThrCorr);
            else
                inds = find(all_corrs_in_compl <= -ThrCorr);
            end
            count = count + max(pos,neg);

            flag = flag || (count < dir_);
            if (count < dir_)
                t(inddir(i1)) = count;
                g(inddir(i1)) = length(unique([matrixGenes(inds);matrixGenes(inddir(i1))])); 
            end
        end
    end
    [dirG,inddirG]=max(g);
    inddir_=find(t==max(t));
    inddirG_=find(g==max(g));
    
    desideredRow=intersect(inddir_,inddirG_);
    if isempty(desideredRow)
        desideredRow=inddirG_;
    end
    if length(desideredRow)==1
        desideredRowFinal=desideredRow;
    end
    if length(desideredRow)>1
        for pp=1:length(desideredRow)
            %colIndex=find(abs(CorrProf(desideredRow(pp),:))>ThrCorr);
            colIndex=find((CorrProf(desideredRow(pp),:))>ThrCorr);
            %tempVect(pp)=sum(abs(CorrProf(desideredRow(pp),colIndex)));
            tempVect(pp)=sum((CorrProf(desideredRow(pp),colIndex)));
            indUnits=find(tempVect==1); %it can be that the same set of profiles have been extracted more than one time for the same complex; the correlation analysis could
            % show some "1". This means
        end
        desideredRowFinal=desideredRow(find(tempVect==max(tempVect)));
        if length(desideredRowFinal)>1
            desideredRowFinal;
            desideredRowFinal=desideredRowFinal(1);
            %here: if desideredRowFinal has more than one value it is
            %necessary to check if it is always the same value of
            %correlation (that means that the correlation is computed 2
            %times on the same couple of profiles...the matrix is
            %simmetric)
        end
%         desideredRowFinalSize=size(desideredRowFinal)
%         if desideredRowFinalSize(2)>1
%             problem=0;
%             clear problem
%         end
%         clear desideredRowFinalSize
    end
   
    
    if dir_ >= 1 
%         [maxt, indt]=max(t);
%         maxdir=maxt;
%         indD=indt;
        indD=desideredRowFinal; % this has to be of size:1x1
        
        %similarRows=find(abs(CorrProf(indD,:))>ThrCorr);
        similarRows=find((CorrProf(indD,:))>ThrCorr);
        if size(similarRows,1)>1
            similarRows=similarRows';
        end
        indprofilesProva=[indD similarRows];
        indToSwap_t=find(CorrProf(indD,:)<-ThrCorr);
        matrixProfiles2(indToSwap_t,:)=matrixProfiles2(indToSwap_t,:)*-1;

        clear maxt indt 

        selectedProfiles=matrixProfiles2(indprofilesProva,:);
        referenceProfile=median(matrixProfiles2(indprofilesProva,:));
        if (simil_type==0)
            tempBestSize=corr([referenceProfile;selectedProfiles]');
        else
            tempBestSize=corr_dist([referenceProfile;selectedProfiles]',0);
        end
        bestSize=length(find(tempBestSize(1,2:end)>ThrCorr));
        clear tempBestSize
        %direction=max(t(j),k(j));
%         fprintf(fid,'#############################\n');
%         fprintf(fid,'Complex \t %i\t NumProteins: %i\n',i,length(newComplexes{i}));
%         fprintf(fid,'Number of total profiles in the complex: %i\n',length(CorrProf));
        ll=length(unique(matrixGenes(indprofilesProva)));
%         percentage=ll/length(newComplexes{i})*100;
%         fprintf(fid,'profile query (row in the corr matrix): %i \t numGenes selected %i\t percentageGenes %2.3f\n',t(indD),ll,percentage);
%         fprintf(fid,'#############################\n\n');
    else
        selectedProfiles=[];
        referenceProfile=[];
        bestSize=0;
        ll=0;
    end
else
    selectedProfiles=[];
    referenceProfile=[];
    bestSize=0;
    ll=0;
end
numSelGenes=ll; clear ll
clear dir_ inddir
clear NprofGene matrixGenes  CorrProf matrixProfiles k t maxdir indD indprofilesProva ll percentage
% matrixProfiles=[];
% matrixGenes=[];
