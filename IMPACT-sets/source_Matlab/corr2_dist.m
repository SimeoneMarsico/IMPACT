function out = corr2_dist(a,b, do_abs)
out = ones(size(a,2),size(b,2));
for i =1:size(a,2)
    for j=1:size(b,2)
        n = norm(a(:,i)-b(:,j));
        if (do_abs)
            n1 = norm(a(:,i)+b(:,j));
            n = min(n,n1);
        end
        if (n == 0)
            n = eps;
        end
        out(i,j) = 1/n;
        
    end
end
