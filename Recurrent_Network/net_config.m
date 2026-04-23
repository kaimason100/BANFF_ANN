clear all
close all
clc
%% Run the initial network configuration for learning. 
rng(1)
%% 
N = 16000; %number of neurons. 
k = 3; %maximum dimension of the dynamical system(s) under question.   
nx = 1000;  %number of sample points for optimization. 
x = (2*rand(k,nx)-1);   %random sample points across space for training. 
xv =  2*rand(k,100)-1; %random sample points across space for validation.
c = 2*rand(1,nx)-1;
cv = 2*rand(1,100)-1;
win = 5*(2*rand(N,1)-1);  %input weight matrix.

NE = round(N/2);
eta = round(2*rand(N,k)-1)*4;
eta = abs(eta);
phi = round((randn(N,k)))/sqrt(N);
phi1(1:NE,:) = abs(phi(1:NE,:));
phi1(NE+1:N,:) = -abs(phi(NE+1:N,:));
phi = phi1;
omega = eta*phi';

tab = table2cell(readtable('dynamics_list.xlsx'));
%%
dynamics = tab(:,1);
dim  = cell2mat(tab(:,2));
input = cell2mat(tab(:,3));
ic = cell2mat(tab(:,4:7));
nsys = length(dynamics);
rescale = cell2mat(tab(:,8));
clear tab
save network_configuration.mat -v7.3  

