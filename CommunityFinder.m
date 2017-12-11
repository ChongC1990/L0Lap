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