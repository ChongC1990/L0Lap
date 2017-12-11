function [ct,out]=hatL0_Lap(A,neta,meta,lambda2)
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

%------------------------------iteration-----------------------------------
while (sum(sum(A))>5)
%     if (length(out.left)<=50)
%         pvalue = CommunityTest(A,A1,neta,meta,allindex,repeat);
%         if(pvalue>0.05)
%             ct{K} = out.left;
%             break
%         end
%     end
    S = hatCommunityFinder(A,A1,neta,meta,lambda2);
    A(S,:) = [];
    A(:,S) = [];
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

function [S] = hatCommunityFinder(A,A1,neta,meta,lambda2)

% Extract one community in pL0Lap

[n,~] = size(A);
[n0,~] = size(A1);
Q = lap(A);
lambda1 = 1/sqrt(n);
H = Q+lambda1*eye(n);
p0 = sum(sum(A1))/((n0^2-n0));
% initialize u,v for inner iteration
v = ones(n,1);
v = v/norm(v,2);
u = rand(n,1);
u=u/norm(u,2);
phi = zeros(1,neta);
    
% value of eta
eta = linspace(0,meta/n,neta);
% inner iteration
S = hatOpt(v,u,H,A,eta(1),lambda2);
phi(1)=Des(S,A,p0);
for i=2:neta
    S = hatOpt(v,u,H,A,eta(1,i),lambda2);
    phi(1,i)=Des(S,A,p0);
%     if (phi(i)<phi(i-1))
%         break
%     end
    %S = indexjust(S, allindex);
    %phi(1,i)=Des(S,A);
    %phi(1,i) = Des(S,A1);
end
    % find the max eta
tmp = find(phi==max(phi));
tmp = tmp(1,1);
%tmp = select(phi);
% tmp = i-1;
[S,u] = hatOpt(v,u,H,A,eta(1,tmp),lambda2);
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

function [S,u]=hatOpt(v,u,H,A,eta,lambda1)
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

while (norm(up-u,2)>1e-4 && k<100)
    up=u;
    u=Alt(H*v,eta);
    vp=v;
    v=Alt(H*u,eta);
    k=k+1;
end
up = zeros(size(u,1),1);
vp = zeros(size(v,1),1);
% Stopping criteria: 
% close enough for u,v between two steps or reaches the max iteration
while (norm(up-u,2)>1e-4 && k<100)
    up = u;
    us = u2us(u,A);
    u=Alt(H*v+lambda1*us,eta);
    vp=v;
    vs = u2us(v,A);
    v=Alt(H*u+lambda1*vs,eta);
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


function u=Alt(z,eta)
% z=Alt(v,H,eta) solves the subproblem: 
%       max{z|z'Hv-eta*norm(z,0)}
% Input:
% v: a vector with unit L2 norm
% H: Graphic Laplician
% eta: Restrict the size of community
% Output:
% z: the result vector with unit L2 norm of subproblem
[n,~]=size(z);
s=zeros(n,1);
u = ones(n,1);
z1=sort(abs(z),'descend');
s(1,1) = norm(z1(1:1),2)-eta*1;
for i=2:n
    s(i,1)=norm(z1(1:i),2)-eta*i;
    % If the previous summation exceeds the next, stop
    if (s(i,1)<s(i-1,1))
        break
    end
end

if i<n
    m=find(abs(z)<z1(i-1,1));
    z(m,1)=0;
    u(m,1) = 0;
end
% the optimal z with unit L2 norm
u=z/norm(z,2);
%u(m,1) = 0;
%u = u/norm(u,2);
end
 
function [us] = u2us(u,A)

% Convert u to us.

[n,~] = size(A);
us = zeros(n,1);
l = find(u~=0);
vol = sum(sum(A(l,:)));
d = sum(A(l,:),2);
us(l) = sqrt(d)/sqrt(vol);
end 
