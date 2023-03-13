function [b,isterminal,direction] = eventfc_line(t,y)
% event function for matlab ode solver which stops integration when
% plume width becomes very large
% this saves time when the plume approaches max height

% plumewidth = y(1)^2/y(2);
% b = sqrt(plumewidth) - 500;  % halt integration when width = 500

velocity = y(2)/y(1);
b = velocity-0.01;

isterminal = 1;        % halt integration 
direction = 0;         % the zero can be approached from either direction