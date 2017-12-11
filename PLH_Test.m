%%%%%%%%%  Pseudo-likelihood Bickel %%%%%%%%%%%%%%%%%%%%%%%%

% PLH
% Model: SBM & DCBM
% without or with outlier
% beta = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2];

%%%%%%%%% without outliers %%%%%%%%%%%%%%%%%%%%%

% SBM
N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data/sbm/PLH.xls',PLH');


% DCBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data/dcbm/PLH.xls',PLH');

%%%%%%%%%%%%%%%%%%%%% with outliers %%%%%%%

% SBM
N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,100,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data1/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data1/sbm/PLH.xls',PLH');

% DCBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,100,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data1/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data1/dcbm/PLH.xls',PLH');

%%%%%%%%%%%%%%%%%%%%%%%%% Changing Lambda %%%%%%%%%%%%%%%%

%SBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data2/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data2/sbm/PLH.xls',PLH');
% DCBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data2/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data2/dcbm/PLH.xls',PLH');


%%%%%%%%%%%%%%%%%%%%%%%%% Changing Lambda %%%%%%%%%%%%%%%%

%SBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:4
    for j=1:N
        name = ['data3/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data3/sbm/PLH.xls',PLH');
% DCBM

N = 1;
net = [ones(100,1); repmat(2,100,1); repmat(3,100,1); repmat(4,100,1); repmat(5,100,1); repmat(6,50,1); repmat(7,50,1);...
    repmat(8,50,1); repmat(9,50,1); repmat(10,50,1); repmat(11,50,1); repmat(12,20,1); repmat(13,20,1); repmat(14,20,1);...
    repmat(15,20,1); repmat(16,20,1); repmat(17,20,1); repmat(18,20,1); repmat(19,20,1); repmat(20,20,1); repmat(21,20,1)];

PLH  = zeros(10,N);
kset = 1:1:25;

for i=1:4
    for j=1:N
        name = ['data3/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k] = community_number_via_spectral_methods(A,kset);
        e = SC(A,k,0.2);
        PLH(i,j) = MuIn(e,net);
    end
end

PLH = mean(PLH,2);    
xlswrite('data3/dcbm/PLH.xls',PLH');