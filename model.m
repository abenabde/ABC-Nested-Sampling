function [t,x]= model(param) 
% The single-degree-of-freedom (SDOF) system governed by the 
% following ordinary differential  equation is considered: 
%                    mx''+cx'+kx=f(t)
global m c k % mass, damping and stiffness parameters
%----------------------------------------------------------
% Input : 
t_start = 0; % starting time
t_end   = 10;  %final time in seconds. 
time_span =[t_start t_end]; 
time_span =t_start:0.01:t_end; 
c=param(1,1);
k=param(1,2);
m = 1; % the mass of the system is supposed to be known 
%----------------------------------------------------------
%Output :
% t : time vector going from t0 to tmax
% x : response matrix [displacement  velocity]
%----------------------------------------------------------
initial_position = 0; 
initial_speed    = 0; 
x0 = [initial_position  initial_speed];  
% Simulate
[t,x]=ode45(@fun,time_span,x0);
end