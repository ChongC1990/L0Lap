function [Sorigin, allindex] = indexjust(S, allindex)
% [Sorigin, allindex] = indexjust(S, allindex) shows the relation between
% the new indexes of S in the sub-matrix of A and the corresponding indexes
% in original A with S.

% Input:
% S: set of index in the sub-matrix of A
% allindex: the remaining indexes in original A
% Output:
% Sorigion: the corresponding indexes in original A with S
% allindex: the new remaining indexes in original A
[n,~]=size(S);
Sorigin = zeros(n,1);
for i=1:n
    Sorigin(i,1) = allindex(S(i,1),1);
end
allindex(S,:)=[];    
end