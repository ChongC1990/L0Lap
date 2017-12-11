function pvalue = cfuse(ct0,A)

%Comupute the p-value for each community in ct0

[~,n] = size(A);
p0 = sum(sum(A))/n/(n-1);
pvalue = zeros(1,length(ct0));

for i=1:length(ct0)
    n0=length(ct0{i});
    X = sum(sum(A(ct0{i},ct0{i})))/2;
    N = n0*(n0-1)/2;
    pvalue(i) = -log10(binocdf(X-1,N,p0,'upper'));
end

end