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