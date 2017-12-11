function [Dense] = confuse(Find,True)

Dense = zeros(max(Find),max(True));

True1=[];
for i=1:max(Find)
    dense = zeros(1,max(True));
    for j=1:max(True)
        dense(j)=length(intersect(find(Find==i),find(True==j)))/length(union(find(Find==i),find(True==j)));
    end
    index = find(dense==max(dense));
    True1 = [True1 index(1)];
    for j=1:max(Find)
        Dense(j,i)=length(intersect(find(Find==j),find(True==index(1))))/length(union(find(Find==j),find(True==index(1))));
    end
end

True1 = setdiff(1:max(True),True1);
for i=1:length(True1)
    for j=1:max(Find)
        Dense(j,max(Find)+i)=length(intersect(find(Find==j),find(True==True1(i))))/length(union(find(Find==j),find(True==True1(i))));
    end
end
end