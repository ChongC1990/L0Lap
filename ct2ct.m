function hct = ct2ct(A,ct,thes,repeat,p0)

% Filter the false community by permutation test

[test,cp] = pfuse(ct,A,thes,repeat,p0);


cp = sort(cp);
hp = cfuse(ct,A);
hct = {};
k=1;
for j=1:length(ct)
    if(hp(j)>=cp(repeat*0.95))
        hct{k}=ct{j};
        k=k+1;
    end
end

if isempty(test)
    hct=ct;
    return
end

end