function label = assign(hct,ct)

label = zeros(1,length(hct));
isc = zeros(1,length(ct));
for i=1:length(hct)
    for j=1:length(ct)
        isc(j) = length(intersect(hct{i},ct{j}));
    end
    label(i)=find(isc==max(isc));
    isc = zeros(1,length(ct));
end