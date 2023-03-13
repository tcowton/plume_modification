%% calculate_plume_outputs
% Tom Cowton 03/23

% Generate a profile of temperature and salinity values of plume modified
% water, in response to varying subglacial discharge

% Uses function gsw_sigma0.m available as part of the Gibbs-SeaWater (GSW) Oceanographic Toolbox
% (https://www.teos-10.org/software.htm)
% McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5. 

% Inputs:
% Q_glacier - subglacial discharge timeseries (m^3/s)
% T_ambient - ambient temperature profile (conservative Temperature, deg C)
% S_ambient - ambient salinity profile (absolute salinity, g/kg)
% z_ambient - depth values for ambient profile
% d_gl - grounding line depth
% p - struture of plume model parameters

% Outputs:
% T_plume - vector of temperature values for plume modified water
% S_plume - vector of salinity values for plume modified water
% z_plume - vector of depth values for plume modified water
% plume_ouputs - matrix containing values output from the plume model for
% each discharge value. Columns refer to:

% 1) Final plume depth (m) (not the same as depth of neutral buoyancy)
% 2) Plume Thickness (m)
% 3) Vertical velocity (m/s)
% 4) Temperature (degC)
% 5) Salinity (g/kg)
% 6) volume flux (m^3/s)
% 7) momentum flux
% 8) melt rate (m/d w.eq.)
% 9) density (kg/m^3)

function[z_plume,T_plume,S_plume,plume_outputs] = calculate_plume_outputs(Q_glacier,T_ambient,S_ambient,z_ambient,d_gl,p)

rho_ambient0 = gsw_sigma0(S_ambient,T_ambient); % ambient potential density anomaly profile

count = 1;
% generate plume outputs
for j = 1:length(Q_glacier)
    if Q_glacier(j)>0 % if there is runoff at that timestep
        
        [A, X, zp] = run_plume_model(Q_glacier(j),T_ambient,S_ambient,-z_ambient,0,-(abs(d_gl)),p); % run plume model
        
        plume_outputs(count,1:9) = [zp(end) X(end,:) A(end,:)];
                
        rho_plume0 = gsw_sigma0(plume_outputs(count,5),plume_outputs(count,4)); % calculate potential density anomaly
        rhodiff = abs(rho_ambient0-rho_plume0);
        plume_outputs(count,10) = z_ambient(find(rhodiff==min(rhodiff),1)); % get depth of neutral buoyancy
        
        count = count+1;
    end
end

% interpolate plume T and S onto regular 1m z spacing based on neutral
% buoyancy depth

[uv,ui] = unique(plume_outputs(:,10)); % get unique values
z_plume = min(plume_outputs(:,10)):max(plume_outputs(:,10));
T_plume = interp1(plume_outputs(ui,10),plume_outputs(ui,4),z_plume);
S_plume = interp1(plume_outputs(ui,10),plume_outputs(ui,5),z_plume);
