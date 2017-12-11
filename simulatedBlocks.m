function [b] = simulatedBlocks(n, net, p1, p2)

m=length(net);
b = zeros(n,n);
P = zeros(n,n);
P(:,:)=p2;
for i=1:m
    k = net{i};
    P(k, k) = p1(i);
end

P = P- diag(diag(P));

for i=1:n
    for j=i:n
        b(i,j)=floor(rand(1) + P(i,j));
    end
end

b = b - diag(diag(b));

for (i=1:n)
    for (j=1:i-1)
        b(i, j) = b(j, i);
    end
end


b = sparse(b);