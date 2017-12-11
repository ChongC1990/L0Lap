function [c,ceq]=confun(x)
c=[];
[n,~]=size(x);
ceq = 0;
for i=1:n
    ceq=ceq+x(i)^2;
end
ceq=ceq-1;
end