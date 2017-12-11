## L0Lap
Community detection by L0-penalized graph laplacian
## Quick Start
``` matlab
net = {1:100, 101:200, 201:300, 301:400, 401:500, 501:550, ...
    551:600, 601:650, 651:700, 701:750, 751:800, 801:820,...
    821:840, 841:860, 861:880, 881:900, 901:920, 921:940, 941:960, 961:980, 981:1000};
n = 1000;
signal = 0.5;
noise = 0.05
A = simulatedBlocks(n,net,signal*ones(1,21),noise);
ct = L0_Lap(A,20,2);
MuIn(net,ct)

