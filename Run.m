%% Parameter estimation using ABC-NS algorithm: SDOF system
clc;clear;close all;
load('training_data.mat');
% true_value = [20 100];
tol = 1000; % initial tolerance threshold 
accuracy=0.01; % change the value of accuacy as desired
theta_ABC = abc_sdof(data,tol,accuracy) 
