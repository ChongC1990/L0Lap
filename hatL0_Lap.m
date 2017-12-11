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
