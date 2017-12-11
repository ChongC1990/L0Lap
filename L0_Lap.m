function [ct,out]=L0_Lap(A,neta,meta)
% [ct,outline,beta]=L0_Lap(A,c,nl) solves the problem of extracting
% different communities for adjacent matrix A. 

% Input: 
% A: Adjacent Matrix
% c: c/n is the max value of eta
% n1: the number of advisable eta
% Output:
% ct: cell aray of the result
% beta£ºthe optimal eta in each iteration

%--------------------------initialization----------------------------------
t = cputime;
A1 = A;
ct = {};
[nfull,~] = size(A);
% Q=lap(A);
% H=Q+lambda*eye(nfull);
%I=eye(n);
K = 1;
m = nfull;
out.left = [1:nfull];
allindex = (1:nfull)';
[nl,~] = size(A);
%------------------------------iteration-----------------------------------
while (nl>=5)
%     if (length(out.left)<=50)
%         pvalue = CommunityTest(A,A1,neta,meta,allindex,repeat);
%         if(pvalue>0.05)
%             ct{K} = out.left;
%             break
%         end
%     end
    S = CommunityFinder(A,A1,neta,meta);
    A(S,:) = [];
    A(:,S) = [];
    [nl,~] = size(A);
    % remove the set S of index from allindex
    [S, allindex] = indexjust(S, allindex);
    %record S as a new community
    out.left = setdiff(out.left,S);
    ct{K} = S;
    K = K+1; 
    m = m-length(S);
end
out.t = cputime-t;
    
end

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

function [S] = CommunityFinder(A,A1,neta,meta)

% Extract one community in L0Lap


[n,~] = size(A);
[n0,~] = size(A1);
Q = lap(A);
lambda = 1/sqrt(n);
H = Q+lambda*eye(n);
p0 = sum(sum(A1))/((n0^2-n0));
% initialize u,v for inner iteration
v = ones(n,1);
v = v/norm(v,2);
u = rand(n,1);
u=u/norm(u,2);
phi = zeros(1,neta);
    
% value of eta
%eta = exp(linspace(-5,0,neta))/n*meta;
eta = linspace(0,meta/n,neta);
% inner iteration
% S = Opt(v,u,H,eta(1));
% phi(1)=Des(S,A,p0);
for i=1:neta
    S = Opt(v,u,H,eta(1,i));
    %S = indexjust(S, allindex);
    phi(1,i)=Des(S,A,p0);
%     if (phi(i)<phi(i-1))
%         break
%     end
    %phi(1,i) = Des(S,A1,A);
end
    % find the max eta
tmp = find(phi==max(phi));
tmp = tmp(1,1);
% tmp = select(phi);
% tmp = i-1;
[S,u] = Opt(v,u,H,eta(1,tmp));
end

function Q=lap(A)
% Q=lap(A) computes the graphic Laplacian Q of adjacent matrix A.

% Input:
% A: n by n symmetric 0-1 adjacent matrix
% Output:
% Q: Graphic Laplacian of A
[n,~]=size(A);
u=sum(A);
u=sqrt(u);
for i=1:n
    u(i)=max(u(i),1);
end
D=diag(u);
Q=inv(D)*A*inv(D);

end

function [S,u]=Opt(v,u,H,eta)
% [S,u]=Opt(v,u,H,eta) solves the main optimization:
%        max{u|u'Qu-eta*norm(u,0)}
% It uses the method of alternating iteration for subproblem:
%        max{u|u'Hv-eta*norm(u,0)-eta*norm(v,0)}
% to get the result.

% Input:
% v: initialize vector v
% u; initialize vector u
% H: H = Q+lambda*eye(n), where Q is the graphic Laplacian of A
% eta: A penalty parameter
% Output:
% S: the index set of community for some eta
% u: the optimal u for the main optimization

k=0;
%while (norm(v-u,2)>1e-4 && k<500)
up = zeros(size(u,1),1);
vp = zeros(size(v,1),1);
% Stopping criteria: 
% close enough for u,v between two steps or reaches the max iteration
while (norm(up-u,2)>1e-4 && k<100)
    up=u;
    u=Alt(H*v,eta);
    vp=v;
    v=Alt(H*u,eta);
    k=k+1;
end
% record the nonzero indexes as S
S=find(u~=0|v~=0);
end

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
