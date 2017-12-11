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
    