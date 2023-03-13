% run_plume_model
% Donald Slater 02/14, updated Tom Cowton 03/23

% Script to solve the buoyant plume model (Jenkins 2011, Slater et al 2016)
% Allows for line plume or point source of buoyancy (half conical plume)

% Uses function gsw_rho.m available as part of the Gibbs-SeaWater (GSW) Oceanographic Toolbox
% (https://www.teos-10.org/software.htm)
% McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5. 

% Inputs:
% Q - subglacial discharge (m^3/s)
% T_ambient - ambient temperature profile (conservative Temperature, deg C)
% S_ambient - ambient salinity profile (absolute salinity, g/kg)
% z_ambient - ambient depth profile (m)
% top - upper limit of plume (m)
% groundingline - glacier grounding line depth (m)
% p - structure of plume parameters

% Outputs:
% A - depth profiles of plume properties for individual plumes, showing how properties evolve as plume rises (i.e. not
% profiles of final plume properties for multiple subglacial discharge
% values)
% 1) volume flux (m^3/s)
% 2) momentum flux
% 3) melt rate (m/d w.eq.)
% 4) density (kg/m^3)
% X - as for A, but containing
% 1) Plume Thickness (m)
% 2) Vertical velocity (m/s)
% 3) Temperature (degC)
% 4) Salinity (g/kg)
% z - depth values relating to A and X


function[A, X, z] = run_plume_model(Q,T_ambient,S_ambient,z_ambient,top,groundingline,p)

% PHYSICAL PARAMETERS
% unpack from structure for convenience

E_0 = p.E_0; % entrainment coefficient
rho0_plume = p.rho0_plume; % reference density of plume (in accordance with Boussinesq)
rho0_ambient=p.rho0_ambient; % reference density (kg/m^3)
g=p.g; % gravity
cw=p.cw; % water heat capacity
ci=p.ci; % ice heat capacity
L= p.L; % latent heat of ice melt
k = p.k; % water-ice drag coefficient
GT= p.GT; % thermal turbulent transfer coefficient (dimensionless)
GS= p.GS; % salt turbulent transfer coefficient (dimensionless)
lambda1=p.lambda1; % local freezing point equation: coefficient of salinity
lambda2=p.lambda2; % local freezing point equation: constant offset
lambda3=p.lambda3; % local freezing point equation: coefficient of depth
Tice=p.Tice; % ice temperature
w = p.w; % plume width (m)
plume_type = p.plume_type; % point or line

% Initial T,S
T_i=lambda2+lambda3*groundingline; % plume temperature (in situ) = pressure melt point
S_i=0.01; % plume salinity

P_ambient = rho0_ambient*g*abs(z_ambient)*1*10^(-4); % pressure in dbar
rho_ambient = zeros(size(T_ambient));

for i=1:length(T_ambient)
    rho_ambient(i)=gsw_rho(S_ambient(i),T_ambient(i),P_ambient(i)); % convert T & S to density
end

% Values used to derive initial u and d
rho_i = gsw_rho(S_i,T_i,groundingline); % intitial plume density
rho_a = rho_ambient(round(z_ambient)==round(groundingline)); % initial ambient density
gprime = g*(rho_a-rho_i)/rho0_plume; % reduced gravity g'

clearvars z X

if strcmp(plume_type,'point') % point source plume
            
    % calcuate u_i and d_i so that there is a balance of buoyancy and momentum 
    d_i = ((32*E_0*Q.^2)/(5*pi.^2*gprime)).^(1/5); % initial plume radius (m)
    u_i = (2/pi)*(((5*pi.^2*gprime)/(32*E_0)).^(2/5))*Q.^(1/5); % initial plume velocity (m/s)
    
    % CALL DIFFERENTIAL EQUATION SOLVER
    opts = odeset('AbsTol',1e-1,'RelTol',1e-3,'Events',@eventfc_pointsource); % ODE solver options
    [z,X]=ode23(@(Z,x) plume(Z,x,p.E_0,g,rho0_plume,k,GT,GS,cw,ci,lambda1,lambda2,lambda3,Tice,L,rho_ambient,T_ambient,S_ambient,z_ambient),[groundingline,top],[d_i,u_i,T_i,S_i],opts);
    
elseif strcmp(plume_type,'line') % line plume
    
%   Calcuate u_i and d_i so that there is a balance of buoyancy and momentum 
    d_i = ((E_0*Q.^2)/(gprime*w.^2)).^(1/3); % initial plume thickness (m)
    u_i = ((gprime*Q)/(E_0*w)).^(1/3); % initial plume velocity (m/s)    
    
    % CALL DIFFERENTIAL EQUATION SOLVER
    opts = odeset('AbsTol',1e-1,'RelTol',1e-3,'Events',@eventfc_line); % ODE solver options
    [z,X]=ode23(@(Z,x) line_plume(Z,x,rho_ambient,T_ambient,S_ambient,z_ambient,E_0,k,rho0_plume,g,cw,ci,L,GT,GS,lambda1,lambda2,lambda3,Tice),[groundingline,top],[d_i,u_i,T_i,S_i],opts);
    
end

% POST PROCESS

matrixsize=size(X);
iterations=matrixsize(1);
A=zeros(iterations,4);

P_z = rho0_ambient*g*abs(z)*1*10^(-4); % pressure in dbar

%% CALCULATE INTERESTING VARIABLES
if strcmp(plume_type,'point')
    for i=1:iterations
        % VOLUME FLUX
        A(i,1)=pi*X(i,1)*X(i,1)*X(i,2)/2;
        % MELT RATE
        a=abs(X(i,2))*lambda1*sqrt(k)*(cw*GT-ci*GS);
        b=-L*GS*sqrt(k)*abs(X(i,2))+ci*GS*sqrt(k)*abs(X(i,2))*lambda1*X(i,4)-ci*GS*sqrt(k)*abs(X(i,2))*(lambda2+lambda3*z(i)-Tice)-cw*GT*sqrt(k)*abs(X(i,2))*(X(i,3)-lambda2-lambda3*z(i));
        c=L*GS*sqrt(k)*abs(X(i,2))*X(i,4)+ci*sqrt(k)*GS*abs(X(i,2))*X(i,4)*(lambda2+lambda3*z(i)-Tice);
        Sb=(1/(2*a))*(-b-sqrt(b^2-4*a*c));
        Tb=lambda1*Sb+lambda2+lambda3*z(i);
        A(i,3)=86400*sqrt(k)*X(i,2)*GS*(X(i,4)-Sb)/Sb;
        % DENSITY
        A(i,4)=gsw_rho(X(i,4),X(i,3),P_z(i));
        % MOMENTUM FLUX
        A(i,2)=pi*X(i,1)*X(i,1)*X(i,2)*X(i,2)*A(i,4)/2;
    end
    
elseif strcmp(plume_type,'line')
    for i=1:iterations
        % volume flux
        A(i,1)=w*X(i,1)*X(i,2);
        % momentum flux
        A(i,2)=X(i,1)*X(i,2)*X(i,2);
        % melt rate
        a=abs(X(i,2))*lambda1*sqrt(k)*(cw*GT-ci*GS);
        b=-L*GS*sqrt(k)*abs(X(i,2))+ci*GS*sqrt(k)*abs(X(i,2))*lambda1*X(i,4)-ci*GS*sqrt(k)*abs(X(i,2))*(lambda2+lambda3*z(i)-Tice)-cw*GT*sqrt(k)*abs(X(i,2))*(X(i,3)-lambda2-lambda3*z(i));
        c=L*GS*sqrt(k)*abs(X(i,2))*X(i,4)+ci*sqrt(k)*GS*abs(X(i,2))*X(i,4)*(lambda2+lambda3*z(i)-Tice);
        Sb=(1/(2*a))*(-b-sqrt(b^2-4*a*c));
        Tb=lambda1*Sb+lambda2+lambda3*z(i);
        A(i,3)=86400*sqrt(k)*X(i,2)*GS*(X(i,4)-Sb)/Sb;
        % density
        A(i,4)=gsw_rho(X(i,4),X(i,3),P_z(i));
    end
end


