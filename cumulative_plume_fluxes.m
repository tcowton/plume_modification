%% plume fluxes
% Tom Cowton 03/23

% Generate a cumulative plume volume flux depth profile, based on plume
% model forced by time series of subglacial discharge

% Inputs:
% plume_fluxes - time series of plume volume flux (corresponds to plume_outputs(:,10))
% plume_nb_depths - time series of plume neutral buoyancy depth (corresponds to plume_outputs(:,6))
% zbins - vector of bounds for depth bins (e.g. 0,10,20... m) in which to accumulate plume fluxes
% dt - time step between plume flux input values (seconds)

% Outputs:
% vp - total cumulative plume volume flux per depth bin (m^3)

function[vp] = cumulative_plume_fluxes(plume_fluxes,plume_nb_depths,zbins,dt)

vp = zeros(length(zbins)-1,1);

for t=1:length(plume_nb_depths)
    for k = 1:length(zbins)-1
        if plume_nb_depths(t)>=zbins(k) && plume_nb_depths(t)<zbins(k+1) % find plume outputs in that depth bin
            vp(k) = vp(k)+plume_fluxes(t)*dt; % convert from m^3/s to m^3/day
        end
    end
end