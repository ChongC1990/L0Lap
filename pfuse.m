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