clc;
clear;
close all;
addpath(genpath('Data'));
load('Ori_WDC.mat');

Ohsi = Img;
Nway = size(Ohsi);
num = round(Nway(3)*0.3);
p_n3 =  0.1+(0.2-0.1)*rand(Nway(3),1);
sigma_n3 =  0.1+(0.2-0.1)*rand(Nway(3),1);
add_band = randperm(Nway(3),num);
save('WDC_level.mat','sigma_n3','p_n3');
save('WDC_band.mat','add_band');