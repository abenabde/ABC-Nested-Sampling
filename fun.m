%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
% Differential equation of the system
%             mx''+ cx' + kx = f(t) 
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
function xdot=fun(t,x) 
global m c k
xdot_1 = x(2); 
xdot_2 = -(c/m)*x(2) - (k/m)*x(1) + force(t)/m; 
xdot = [xdot_1 ; xdot_2]; 
end 