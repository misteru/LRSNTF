clc;
clear;
close all;
addpath(genpath('data'));
addpath(genpath('Noises'));

i=1;
methodname  = { 'Noise', 'LRSNTF'};
Mnum = length(methodname);

case_num=1;
load(strcat('WDC_case',num2str(case_num),'.mat'))

Nway = size(Ohsi);


Re_hsi  =  cell(Mnum,1);
psnr    =  zeros(Mnum,1);
ssim    =  zeros(Mnum,1);
sam     =  zeros(Mnum,1);
time    =  zeros(Mnum,1);

addpath(genpath('lib'));
addpath(genpath('data'));
opts=[];
opts.R       = 8;
opts.rho     = 0.1;
opts.tau     = 0.2;
opts.lambda  = 1000;
opts.beta    = 1500;
opts.mu      = 0.04;
opts.max_it  = 60;
opts.Bmax_it = 10;
opts.tol     = 1e-4;

fprintf('\n');
t0= tic;
[Re_hsi{i},A,B,S,Out] = LRSNTF(Nhsi, opts);
time(i) = toc(t0);
[psnr(i), ssim(i), sam(i)] = HSIQA(Ohsi * 255, Re_hsi{i} * 255);
fprintf('Ours: PSNR = %5.4f  time = %5.2f\n',  psnr(i),time(i));


