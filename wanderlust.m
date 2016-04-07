addpath(genpath('cyt/src'))

data = dlmread('wanderlust/data/HSMM_E.txt');
data = data';
[traj,lnn] = wanderlust(data,'cosine',5,30,25,216,20,true,[]);
traj = mean(traj, 1);
traj = traj';
[sorted, order] = sort(traj);
order = fliplr(order);
fileID = fopen('wanderlust/output/HSMM.txt','w');
fprintf(fileID,'%i\n',order);
fclose(fileID);

data = dlmread('wanderlust/data/HSMMgene_Y.txt');
data = data';
[traj,lnn] = wanderlust(data,'cosine',5,30,25,216,20,true,[]);
traj = mean(traj, 1);
traj = traj';
[sorted, order] = sort(traj);
order = fliplr(order);
fileID = fopen('wanderlust/output/HSMM_gene.txt','w');
fprintf(fileID,'%i\n',order);
fclose(fileID);

data = dlmread('wanderlust/data/LPS_E.txt');
data = data';
[traj,lnn] = wanderlust(data,'cosine',5,30,25,122,20,true,[]);
traj = mean(traj, 1);
traj = traj';
[sorted, order] = sort(traj);
order = fliplr(order);
fileID = fopen('wanderlust/output/LPS.txt','w');
fprintf(fileID,'%i\n',order);
fclose(fileID);

data = dlmread('wanderlust/data/qNSC_E.txt');
data = data';
[traj,lnn] = wanderlust(data,'cosine',5,30,25,55,20,true,[]);
traj = mean(traj, 1);
traj = traj';
[sorted, order] = sort(traj);
order = fliplr(order);
fileID = fopen('wanderlust/output/qNSC.txt','w');
fprintf(fileID,'%i\n',order);
fclose(fileID);