% Load in data

function A=LoadData(str)
% str: name of data
% A: Adjacent Matrix

data=load(str);
n=max(max(data));
A=sparse(data(:, 1), data(:, 2), 1, n, n);
A = max(A, A');

end