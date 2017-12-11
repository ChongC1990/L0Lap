function cm = confusionMatrix (e1, e2)

% cm = confusionMatrix (clust1, clust2)
% compute the confusion matrix of two clustering results

cm = zeros(max(e1), max(e2));

for i = 1:max(e1)
    for j = 1:max(e2)
        cm(i, j) = length(intersect(find(e1==i), find(e2==j)));
    end
end
