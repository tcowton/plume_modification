% Set plume parameter values

p.rho0_plume=1000; % reference density of plume (in accordance with Boussinesq)
p.rho0_ambient=1027; % reference density of ambient seawater
p.g=9.81; % gravity
p.cw=3974; % water heat capacity
p.ci=2009; % ice heat capacity
p.L=334000; % latent heat of ice melt

p.k=0.0025; % water-ice drag coefficient
p.GT=0.022; % thermal turbulent transfer coefficient (dimensionless)
p.GS=0.00062; % salt turbulent transfer coefficient (dimensionless)

p.lambda1=-0.0573; % local freezing point equation: coefficient of salinity
p.lambda2=0.0832; % local freezing point equation: constant offset
p.lambda3=0.000761; % local freezing point equation: coefficient of depth
p.Tice=-10; % ice temperature