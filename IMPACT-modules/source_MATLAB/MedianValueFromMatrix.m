function OutValue=MedianValueFromMatrix(CorrMatrix)
OutValue=zeros(1,5);


if size(CorrMatrix,1)>1
    [indCorr_values]=find(CorrMatrix~=0);
    median_r2temp=median(CorrMatrix(indCorr_values).^2); %median of R2 values
    OutValue(1)=median_r2temp;

else
    OutValue(1)=0;
end


if size(CorrMatrix,1)>1
    corr_vector=CorrMatrix(1,2:end);
    [indCorr_values]=find(corr_vector~=0);
    tempCorr=corr_vector(indCorr_values);
    tempCorr=corr_vector(corr_vector~=0);
    mean_r2temp=mean(tempCorr.^2); %mean of R2 values
    OutValue(2)=mean_r2temp;

else
    OutValue(2)=0;
end

if size(CorrMatrix,1)>1
    corr_vector=CorrMatrix(1,2:end);
    mean_temp=mean(abs(corr_vector(corr_vector~=0)));
    OutValue(3)=mean_temp;

else
    OutValue(3)=0;
end

if size(CorrMatrix,1)>1
    corr_vector=CorrMatrix(1,2:end);
    sumR_temp=sum(abs(corr_vector(corr_vector~=0)));
    OutValue(4)=sumR_temp;

else
    OutValue(4)=0;
end

if size(CorrMatrix,1)>1
    corr_vector=CorrMatrix(1,2:end);
    sumR2_temp=sum(corr_vector(corr_vector~=0).^2);
    OutValue(5)=sumR2_temp;

else
    OutValue(5)=0;
end


