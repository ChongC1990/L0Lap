function mi = MuIn(c1,c2)

%Compute the normalized mutual information
%c1 and c2 are two community structure

cm = confusionMatrix(c1,c2);
n = sum(sum(cm));
cm = cm/n;
p1 = sum(cm,2);
p2 = sum(cm,1);

I = 0;
H1 = 0;
H2 = 0;
for i=1:length(c1)
    for j=1:length(c2)
        if (cm(i,j)~=0)
            I = I+cm(i,j)*log(cm(i,j)/p1(i)/p2(j));
        end
    end
end

for i=1:length(c1)
    if(p1(i)~=0)
        H1 = H1-p1(i)*log(p1(i));
    end
end

for i=1:length(c2)
    if(p2(i)~=0)
        H2 = H2-p2(i)*log(p2(i));
    end
end

mi = 2*I/(H1+H2);
end


function cm = confusionMatrix (e1, e2)

% cm = confusionMatrix (clust1, clust2)
% compute the confusion matrix of two clustering results

cm = zeros(length(e1), length(e2));

for i = 1:length(e1)
    for j = 1:length(e2)
        cm(i, j) = length(intersect(e1{i}, e2{j}));
    end
end
end
