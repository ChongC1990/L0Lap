function [us] = u2us(u,A)

% Convert u to us.

[n,~] = size(A);
us = zeros(n,1);
l = find(u~=0);
vol = sum(sum(A(l,:)));
d = sum(A(l,:),2);
us(l) = sqrt(d)/sqrt(vol);
end

