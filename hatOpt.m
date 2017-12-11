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