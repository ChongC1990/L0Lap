function hct = PermutationFilter(A,ct,thes,repeat,p0)

% Filter the false community by permutation test

[test,cp] = pfuse(ct,A,thes,repeat,p0);


cp = sort(cp);
hp = cfuse(ct,A);
hct = {};
k=1;
for j=1:length(ct)
    if(hp(j)>=cp(repeat*0.95))
        hct{k}=ct{j};
        k=k+1;
    end
end

if isempty(test)
    hct=ct;
    return
end

end

function [test,pvalue] = pfuse(ct0,A,thes,repeat,p0)

% Generete reference p-value in permutation test

test = [];
pvalue = ones(1,repeat);
for i=1:length(ct0)
     if(length(ct0{i})<=thes)
         test =[test;ct0{i}];
     end
end
if isempty(test)
    return;
end
pvalue = sizetest(A(test,test),repeat,p0);
end

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
