function pvalue = sizetest(A,repeat,p1)

% Generete reference p-value in permutation test

[nnode,~] = size(A);
nedge = sum(sum(A))/2;
pvalue = zeros(1,repeat);
p0 = (2*nedge)/(nnode^2-nnode);

for i=1:repeat
    A0 = gmatrix(nnode,nedge);
    S = CommunityFinder(A0,A0,10,1);
    n0=length(S);
    X = sum(sum(A0(S,S)))/2;
    N = n0*(n0-1)/2;
    pvalue(i) = -log10(binocdf(X-1,N,p0,'upper'));
end

end