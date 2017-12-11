function showClusters(clusts, a, j, noLine, color, width)
if (nargin < 4)
    noLine = 0;
end

if (nargin < 5)
    color = 'r';
end

if (nargin < 6)
    width = .5;
end

figure (j);
o=[];
for i=1:length(clusts)
%    c = clusts{i};
%    c = c(randperm(length(c)));
   o = [o; clusts{i}];  
%o = [o c];
end
imagesc(1-a(o, o));
colormap(gray);
l = 0.5;
n=length(o) + .5;
       line([0, n], [l, l], 'Color', color, 'LineWidth', width);
       line([l, l], [0, n], 'Color', color, 'LineWidth', width);   
for i=1:(length(clusts))
   l = l + length(clusts{i});
   if (~noLine)
       line([0, n], [l, l], 'Color', color, 'LineWidth', width);
       line([l, l], [0, n], 'Color', color, 'LineWidth', width);   
%       line([0, n], [l, l], 'LineWidth', 1);
%       line([l, l], [0, n], 'LineWidth', 1);   
   end
end
