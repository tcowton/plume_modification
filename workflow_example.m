%% workflow_example

% This script runs through the functions required to calculate the plume
% fluxes and water fractions for a given fjord and year based CTD and
% runoff inputs

% It requires the scripts gsw_sigma0.m and gsw_rho.m available as part of
% the Gibbs-SeaWater (GSW) Oceanographic Toolbox
% (https://www.teos-10.org/software.htm) (McDougall, T.J. and P.M. Barker,
% 2011)

% It includes example CTD data from Kangerlussuup Sermia in 2019, obtained
% from NASA's Ocean's Melting Greenland project archive (OMG 2020)
% (https://podaac.jpl.nasa.gov/omg), and timeseries of subglacial discharge
% based on daily ice sheet surface runoff outputs from the regional climate
% model RACMO2.3, statistically downscaled to 1 km resolution (Noël et al.,
% 2016)

% References:

% McDougall, T.J. and P.M. Barker,
% 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW)
% Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.

% Noël, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I.,
% Fettweis, X., & Van Den Broeke, M. R. (2016). A daily, 1 km resolution
% data set of downscaled Greenland ice sheet surface mass balance
% (1958–2015). The Cryosphere, 10(5), 2361-2377.

% OMG. (2020). OMG CTD Conductivity Temperature Depth. Ver. 1. PO.DAAC, CA,
% USA. https://doi.org/https://doi.org/10.5067/OMGEV-CTDS1

%%  add paths
% add path to Gibbs-SeaWater (GSW) Oceanographic Toolbox

addpath(genpath('C:\Users\trc6\OneDrive - University of St Andrews\Research\Matlab scripts\gsw_matlab_v3_06_13'))

%% load data

load CTD_KS_2019.mat % pre-processed CTD data

% Structures are:
% T.s, S.s and z.s - conservative temperature (degrees C), absolute salinity (g/kg), and depth
% (m) from shelf CTD data
% T.f, S.f and z.f - conservative temperature (degrees C), absolute salinity (g/kg), and depth
% (m) from fjord CTD data

load Q_KS_2019.mat % pre-processed runoff data

% Q_glacier - runoff at daily intervals, 2019 (m^3/s)
% t_matlab - corresponding time values (matlab date format)

%% Define model parameters

% for plume model
p.plume_type = 'line'; % specify 'point' or 'line' plume geometry
p.E_0 = 0.1; % entrainment coefficient
p.w = 200; % plume width (m) (for line plume only)
plume_parameters % load other parameters

d_gl = 250; % grounding line depth (m)

resn_threshold = 1e-2; % threshold for residuals for mixing model (sum of squares)

zfrac = 5:10:max(z.f); % depth intervals for mixing model
zbin = 0:10:max(z.f); % depth bands for cumulative plume fluxes

%% main calculations

% plume model
[z.p,T.p,S.p,plume_outputs] = calculate_plume_outputs(Q_glacier,T.s,S.s,z.s,d_gl,p);

% calculate cumulative plume volume flux
vp = cumulative_plume_fluxes(plume_outputs(:,6),plume_outputs(:,10),zbin,86400);

% calculate potential densities
rho0.f = gsw_sigma0(S.f,T.f);
rho0.s = gsw_sigma0(S.s,T.s);
rho0.p = gsw_sigma0(S.p,T.p);

% mixing model
[fs,fp,zf,zs,zp,xf,xs,xp,resn,resi] = mixing_model(S.f,T.f,rho0.f,z.f,S.s,T.s,rho0.s,z.s,S.p,T.p,rho0.p,z.p,zfrac);
                                             
% after mixing model, apply filter based on residuals, and smooth data
fpf = fp.*(resn<resn_threshold);
fpfs = movmean(fpf,3,'omitnan');

%% plot

ylims = [0 500];

figure
subplot(1,3,1)
hold on; box on; axis ij
plot(T.s,z.s)
plot(T.f,z.f)
plot(T.p,z.p)
xlabel('Temperature (\circC)')
ylabel('Depth (m)')
legend('Shelf','Fjord','Plume','Location','southeast')
ylim(ylims)

subplot(1,3,2)
hold on; box on; axis ij
barh(zbin(1:end-1)+diff(zbin(1:2)),vp/1000.^3)
xlabel('Volume flux (km^3)')
ylabel('Depth (m)')
ylim(ylims)

subplot(1,3,3)
hold on; box on; axis ij
plot(fpfs,zfrac)
xlabel('Fraction')
ylabel('Depth (m)')
ylim(ylims)
