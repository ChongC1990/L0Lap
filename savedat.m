function savedat(name,A)

f=fopen(name,'a');
[~,n]=size(A);
for i=1:(n-1)
    for j=(i+1):n
        if (A(i,j)==1)
            fprintf(f,'%d\t',i);
            fprintf(f,'%d\n',j);
        end
    end
end
fclose(f);
end