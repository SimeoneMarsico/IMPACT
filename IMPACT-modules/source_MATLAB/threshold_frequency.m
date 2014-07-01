function freq = threshold_frequency (AllPROFnet,SeedNodePatternOriginal,n,SIMIL_THR,n_min, N_max, simil_type)
% number of reference profiles is always 1
% number of genes for each group n 

 
N_p = zeros(size(AllPROFnet));
for i = 1 : length(AllPROFnet)
    N_p(i) = size(AllPROFnet{i},1);
end
if (length(find(N_p == 1)) == length(N_p))  % 1 profile per gene, as in the case of average/median/mode profile
    N_max = 1;
end
% --------------------------------------------------------------
% automated assessment of bin intervals
temp = sort(N_p);
for i = 1: N_max
    N(i) = temp(i * floor(length(temp)   /(N_max+1)));
end
N = unique(N);
N_max = min(N_max, length(N));

% binning the dataset in bins
Index = {};
Choise_g = {};
for i = 1 : N_max
    if i == 1
        Index{i} = find(N_p <= N(i));
    elseif i < N_max
        Index{i} = find(N_p > N(i-1) &  N_p <= N(i));
    else
        Index{i} = find(N_p >= N(i));
    end
    numChoise_g(i) = min(n,length(Index{i}));
    Choise_g{i} = Index{i}(randi_nr(length(Index{i}),numChoise_g(i)));
end
% --------------------------------------------------------------
freq = {};
for j = 1 : N_max
   freq{j} = p_vs_g(SeedNodePatternOriginal,AllPROFnet(Choise_g{j}),SIMIL_THR,n_min,simil_type);
end
% --------------------------------------------------------------
for i = 1 : N_max
    freq{i} = freq{i} ./ numChoise_g(i);
end
 

%% -------------------------------------------------------------
function R = randi_nr (n,m)
R = randperm(n);
R = R(1:m);
%% -------------------------------------------------------------
function S = p_vs_g (p,g,t,n_min,simil_type)
S = 0;
for i = 1 : length(g)
    if (simil_type == 0)
        a_corr = (corr(p',(g{i})'));
    else
        a_corr = corr2_dist(p',(g{i})',0);
    end
    T1 = (a_corr >=  t);
    T2 = (a_corr <= -t);
    if (n_min < 1)
        value = max(sum(sum(T1,2) >= (n_min * size(g{i},1))), sum(sum(T2,2) >= (n_min * size(g{i},1))));
        S = S + value; % check number of profiles above a certain fraction
    else
        value = max(sum(sum(T1,2) >= n_min), sum(sum(T2,2) >= n_min));
        S = S + value; % check number of profiles above n_min
    end
end