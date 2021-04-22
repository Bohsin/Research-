function RelaxationPaths

integer_partitions(10,10)

end







function S = integer_partitions(n,count)
%integer partitions of n
if nargin == 1
    count = n;
end
if n < 0 || n ~= round(n)
    error('Only nonnegative integers allowed!');
elseif n == 0
    if count == 0
        S = 0;
    else
        S = zeros(1,count);
    end
else
    x = ones(1,n);
    x(1) = n;
    m = 1;
    h = 1;
    M = [x(1:m) zeros(1,n-m)];
    while x(1) ~= 1
        if x(h) == 2 
           m = m + 1;
           x(h) = 1;
           h = h - 1;
        else
           r = x(h) - 1;
           t = m - h + 1;
           x(h) = r;
           while t >= r
               h = h + 1;
               x(h) = r;
               t = t - r;
           end
           if t == 0
               m = h;
           else
               m = h + 1;
               if t > 1
                   h = h + 1;
                   x(h) = t;
               end
           end
        end
        M = cat(1,M,[x(1:m) zeros(1,n-m)]);
    end
    if count > n
        M = cat(2,M,zeros(size(M,1),count-n));
    end
    S = [];
    for i = 1:size(M,1)
        if(sum(M(i,1:count)) == n)
            S = cat(1,S,M(i,1:count));
        end
    end
end
end