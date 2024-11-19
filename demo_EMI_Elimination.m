close all
clear all
clc

addpath('.\data');
addpath('.\calib');
addpath('.\func');

%% Data preparation
ksp = [];
for i = 1:5
    load(['sample_part' num2str(i) '.mat']);
    ksp = cat(5, ksp, tmp);
    clear tmp;
end
seq = 2.1;
ksp = double(ksp);
para = SeqPara(seq, 1);
para.seq.Nexp = 1;
para.Npri = 1;

data = ksp;
[nx, ny, nz, nex, nc] = size(data);
data = reshape(data, [nx, ny*nz*nex, nc]);

%% Extract the peripheral k-space for calibration
data_calib = circshift(ksp, [0 round(ny/2) round(nz/2) 0 0]);
My= 26; Mz = 5;
data_calib = crop(data_calib, [nx, ny-My, nz-Mz, nex, nc]);
data_calib = reshape(data_calib, [nx, (ny-My)*(nz-Mz)*nex, nc]);
data = fftc(data, 1); data_calib = fftc(data_calib, 1);

%% EMI elimination
EMI_Elimination_Method = 1;
switch EMI_Elimination_Method
    case 1 % rEDITER
        para.reg = 10^(-8);
        para.flag_spectrum = 0;
        [data_emi, ~] = EMI_EDITER(data_calib, data, para);
    case 2 % sEDITER
        para.reg = 10^(-8);
        para.flag_spectrum = 1;
        [data_emi, ~] = EMI_EDITER(data_calib, data, para);
    case 3 % TEE
        [data_emi, ~] = EMI_TEE(data_calib, data, 1);
    case 4 % Without EMI elimination
        data_emi = data(:, :, 1)*0;
end
data_emi = fftc(data_emi, 1);
data = data(:, :, 1)-data_emi;
data = ifftc(data, 1);
data = reshape(data, [nx, ny, nz, nex]);

%% Combination for all averages
tmp = data;
[nx, ny, nz, nex] = size(tmp);
tmp = reshape(tmp, [nx*ny*nz, nex]);
[u, s, v] = svd(tmp, 'econ');
ind_keep = 1;
tmp = u(:, ind_keep)*s(ind_keep, ind_keep)*v(:, ind_keep)';
data = reshape(tmp(:, 1), [nx, ny, nz]);

%% Image reconstruction
data = KSP_ZF_ellipse(data, para);
data = ifft3c(data);
img = abs(data);
img = img(:, :, 13:2:19);
img = img2montage(img, [1 4], 2);
img = img/prctile(img(:), 99);
imshow(img);


