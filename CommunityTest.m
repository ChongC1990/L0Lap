function [pvalue] = CommunityTest(A,A1,neta,meta,allindex,repeat)

[nnode,~] = size(A);
nedge = sum(sum(A))/2;
De = zeros(1,repeat);
De0 = Des(CommunityFinder(A,A1,neta,meta,allindex),A);
for i=1:repeat
    A0 = gmatrix(nnode,nedge);
    S = CommunityFinder1(A0,neta,meta);
    De(i) = Des(S,A0);
end

pvalue = length(find(De>De0))/repeat;

end