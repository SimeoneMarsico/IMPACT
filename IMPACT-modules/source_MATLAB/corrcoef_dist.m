function out = corrcoef_dist(a, do_abs)
out = ones(size(a,2),size(a,2));
for i =1:size(a,2)
    for j=(i+1):size(a,2)
        n = norm(a(:,i)-a(:,j));
        if (do_abs)
            n1 = norm(a(:,i)+a(:,j));
            n = min(n,n1);
        end
        if (n == 0)
            n = eps;
        end
        out(i,j) = 1/n;
        out(j,i) = out(i,j);
    end
end
