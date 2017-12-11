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