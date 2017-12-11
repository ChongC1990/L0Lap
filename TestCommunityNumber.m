% Method: NB, BH, PLH
% Model: SBM & DCBM
% without or with outlier
% beta = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2];

%%%%%%%%% without outliers %%%%%%%%%%%%%%%%%%%%%

% SBM
N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data/sbm/cmn.xls',cmn);


% DCBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data/dcbm/cmn.xls',cmn);

%%%%%%%%%%%%%%%%%%%%% with outliers %%%%%%%

% SBM
N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data1/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data1/sbm/cmn.xls',cmn);

% DCBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data1/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data1/dcbm/cmn.xls',cmn);

%%%%%%%%%%%%%%%%%%%%%%%%% Changing Lambda %%%%%%%%%%%%%%%%

%SBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data2/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data2/sbm/cmn.xls',cmn);

% DCBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:10
    for j=1:N
        name = ['data2/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,10);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data2/dcbm/cmn.xls',cmn);



%%%%%%%%%%%%%%%%%%%%%%%%% Changing Lambda %%%%%%%%%%%%%%%%

%SBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:4
    for j=1:N
        name = ['data3/sbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,N);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data3/sbm/cmn.xls',cmn);

% DCBM

N = 10;
cnum = {};
cnum{1}{1}=zeros(4,N);
cnum{2}{1}=zeros(4,N);
cnum{3}{1}=zeros(4,N);
kset = 1:1:25;

for i=1:4
    for j=1:N
        name = ['data3/dcbm/','b',num2str(i),'r',num2str(j),'.dat'];
        A = LoadData(name);
        [k1,k2,k3] = community_number_via_spectral_methods(A,kset);
        cnum{1}{1}(i,j) = k1;
        cnum{2}{1}(i,j) = k2;
        cnum{3}{1}(i,j) = k3;
    end
end

cmn=zeros(3,4);
cmn(1,:)=mean(cnum{1}{1},2);    
cmn(2,:)=mean(cnum{2}{1},2);
cmn(3,:)=mean(cnum{3}{1},2);

xlswrite('data3/dcbm/cmn.xls',cmn);

