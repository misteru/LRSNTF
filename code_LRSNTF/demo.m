clc;
clear;
close all;
addpath(genpath('lib'));
addpath(genpath('Noises'));

methodname = {'Noise','LRSNTF'};
Mnum = length(methodname);

case_num=5;
load(strcat('WDC_case',num2str(case_num),'.mat'))

%% evaluation indexes
Re_hsi  =  cell(Mnum,1);
psnr    =  zeros(Mnum,1);
ssim    =  zeros(Mnum,1);
sam     =  zeros(Mnum,1);
time    =  zeros(Mnum,1);

%%
Nway = size(Ohsi);

%% corrupted image
i  = 1;
Re_hsi{i} = Nhsi;
[psnr(i), ssim(i), sam(i)] = HSIQA(Ohsi * 255, Re_hsi{i} * 255);
enList = 1;

%% Performing LRSNTF
i = i + 1;
opts=[];
opts.R       = 8;
opts.rho     = 0.1;
opts.tau     = 0.2;
opts.lambda  = 1000;
opts.beta    = 1500;
opts.mu      = 0.04;
opts.max_it  = 50;%50;
opts.Bmax_it = 10;
opts.tol     = 1e-4;

opts.Xtrue     = Ohsi;

%%%%%
fprintf('\n');
disp(['performing ',methodname{i}, ' ... ']);
t0= tic;
[Re_hsi{i},A,B,S,Out] = LRSNTF(Nhsi, opts);
time(i) = toc(t0);
[psnr(i), ssim(i), sam(i)] = HSIQA(Ohsi * 255, Re_hsi{i} * 255);
fprintf('OurMethod8 : PSNR = %5.4f  time = %5.2f\n',  psnr(i),time(i));
enList = [enList,i];


%% Show result
fprintf('\n');
fprintf('================== Result ==================\n');
fprintf(' %8.8s    %5.4s      %5.4s    %5.4s    \n', 'method','PSNR', 'SSIM', 'SAM');
for i = 1:length(enList)
    fprintf(' %8.8s    %5.4f    %5.4f    %5.4f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), sam(enList(i)));
end
fprintf('================== Result ==================\n');

%%
figure,
showHSIResult(Re_hsi,Ohsi,0,1,methodname,enList,1,Nway(3))
