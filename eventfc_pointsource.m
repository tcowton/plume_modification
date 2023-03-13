function [b,isterminal,direction] = eventfc_pointsource(t,y)
% event function for matlab ode solver which stops integration when
% plume width becomes very large
% this saves time when the plume approaches max height

b = y(2)-0.05; % halt integration when plume velocity drops below this level

isterminal = 1;        % halt integration 
direction = 0;         % the zero can be approached from either direction