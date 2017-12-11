function w=Des(S,A,p0)
% w=Des(S,A) calculates the ratio of pin:psum (the psum may be adjusted for
% different problem), which indicates the tightness of selection S in 
% adjacent matrix A.

%Input:
% S: a subset of index
% A: Adjacent Matrix
% Output:
% w: ratio for tightness
[m,n] = size(A);
S1 = [1:1:m]';
S2 = setdiff(S1,S);
os = sum(sum(A(S,S)));
bs = sum(sum(A(S,S2)));
[s,~]=size(S);
pout = 0;
pin = os/(s*s-s);
if s >=m-1
    pout=p0;
else
    pout=bs/(s*(m-s));
end
w = pin/(pin+pout);
%w = pin-pout;
%w = pin/(pin+pout+1);
end