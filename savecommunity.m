function savecommunity(filename,ctT,ctF,id)

n = length(id);
A = cell(n+1,3);
A(:,1)={'id';id};
l1=length(ctT);
l2=length(ctF);
for i=1:l1
    A(ctT{i}+1,2)={i};
end

for i=1:l2
    A(ctF{i}+1,3)={i};
end
B(1,1)='id';
B(1,2)='T';
B(1,3)='F';
xlswrite(filename,A);

end