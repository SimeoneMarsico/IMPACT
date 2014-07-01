function p = gene_network_random (AllPROFnet,AdjMatrixNet,SeedNodeUsed,SeedPatternUsed,SeedPatternUsedReal,moduleSet,PReal,t, N_max,simil_type)

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



P = cell(size(PReal));
for i = 1 : size(SeedPatternUsed,1)
    if size(SeedPatternUsed)==size(SeedPatternUsedReal)
        P{i} = PReal{i};
    else
        if(simil_type==0)
            CC = corrcoef([SeedPatternUsed(i,:); SeedPatternUsedReal]');
        else
            CC = corrcoef_dist([SeedPatternUsed(i,:); SeedPatternUsedReal]',0);
        end
        [maxCC,indmaxCC] = max(CC(1,2:end));
        P{i} = PReal{indmaxCC(1)};
    end
end
    
p = ones(size(moduleSet));

t10 = t * 10; 

for i = 1 : length(moduleSet)
   
    Seed = SeedNodeUsed(i); % find the seed gene
    Q = cell2mat(P{i}'); P_t = Q(:,1);
                                         
    neighbours = find(AdjMatrixNet(Seed,:) == 1); % the current neighbours
    moduleNodes = intersect(moduleSet{i},neighbours); % the module nodes in
                                                      % the current neighbours
    AllNodes = [Seed neighbours]; % all the nodes including all the neighbours
                                  % till the ongoing step
    j = 1; % 'j' serves as the counter, which notes how many genes in 
           % the module set has been considered
    k = length(moduleNodes);
    n = length(neighbours);
    BC = nchoosek(n,k); % binomial coefficient
    S = neighbourhood(AllPROFnet,neighbours, N, N_max); % S shows the numbers of profiles 
                                              % of the current neighbour genes
    PS_t = P_t(S); % the vector of P's with respect to the current 
                   % neighbourhood
    if (size(PS_t,1) < size(PS_t,2))
        PS_t = PS_t';
    end
    if ~isempty(S) && BC <= 10000 % if the current neighbourhood is not empty
                                  % and the binomial coefficient is less
                                  % equal than 10000 (the threshold), we
                                  % use the accurate computation
        
        pp = 1 - dbinocdf(k-1,n,PS_t);
    elseif ~isempty(S) && BC > 10000 % if the current neighbourhood is not empty
                                     % and the binomial coefficient is
                                     % greater than 10000, we consider the
                                     % approximation
        pp = 1 - binocdf(k-1,n,mean(PS_t));
    else % if the current neighbourhood is empty, we force the probability to 1
        pp = 1;
    end
    p(i) = p(i) * pp;
    while j < length(moduleSet{i}) % if not all of the nodes in the module
                                   % set, the process continues
        neighbours = setdiff(find(sum(AdjMatrixNet(moduleNodes,:),1) ~= 0),AllNodes);
        moduleNodes = intersect(moduleSet{i},neighbours);
        AllNodes = [AllNodes neighbours];
        j = j + k; % the new 'j' is obtain by 'j + k', because k new nodes
                   % in the module set have been considered
        k = length(moduleNodes);
        n = length(neighbours);
        BC = nchoosek(n,k);
        S = neighbourhood(AllPROFnet,neighbours, N, N_max);
        PS_t = P_t(S);
        if ~isempty(S) && BC <= 10000
           pp = 1 - dbinocdf(k-1,n,PS_t);
        elseif ~isempty(S) && BC > 10000
            pp = 1 - binocdf(k-1,n,mean(PS_t));
        else
            pp = 1;
        end
        p(i) = p(i) * pp;
    end
end


function p = pp (P,v,w)
p = prod([P(v); 1-P(w)]);



function p = dbinocdf (k,n,P)
if (size(P,1) < size(P,2))
    P = P';
end
p = 0;
for i = 0 : k
    v = nchoosek(1:n,i);
    for j = 1 : size(v,1)
        p = p + pp(P,v(j,:),setdiff(1:n,v(j,:)));
    end
end


function S = neighbourhood (AllPROFnet,neighbours, N, N_max)

if ~isempty(neighbours)
    for i = 1 : length(neighbours)
        N_p = size(AllPROFnet{neighbours(i)},1);
        ind = 1;
        for i2 = 1 : N_max
            if i2 == 1
                if (~isempty(find(N_p <= N(i2))))
                    ind = i2;
                end
            elseif i2 < N_max
                if (~isempty(find(N_p > N(i2-1) &  N_p <= N(i2))))
                    ind = i2;
                end
            else
                if (~isempty(find(N_p >= N(i2))))
                    ind = i2;
                end
            end

        end
    
        S(i) = ind;
        %S_temp = min(size(AllPROFnet{neighbours(i)},1),9);
        %S(i) = min(max(4,size(AllPROFnet{neighbours(i)},1)),7) - 3;
    end
else
    S = [];
end