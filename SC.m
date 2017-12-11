function [e] = SC(A,k,eps)
%spectral clustering 
if k==1
    e = ones(size(A,2),1);
else
    d = sum(A,2); d = d + eps*mean(d); d = 1./sqrt(d);
    B = diag(d) * A * diag(d);
    [U,~] = eigs(B,k);
    e = kmeans(U,k,'replicates',10,...
        'onlinephase','off','emptyaction','singleton');
end
end