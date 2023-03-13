%% mixing_model
% Tom Cowton, 03/23

% Calculate fractions of shelf water and PMW required to give best fit
% to fjord water properties

% Inputs:
% S_fjord, T_fjord, rho0_fjord, z_fjord - profiles of fjord absolute
% salinity (g/kg), conservative temperature (deg C), potential density
% anomaly and depth
% S_shelf, T_shelf, rho0_shelf, z_shelf - as above, for shelf
% S_plume, T_plume, rho0_plume, z_plume - as above, for plume outputs
% z_fractions - depth intervals at which to calculate fractions (m)

% Outputs:
% fs - fraction shelf water
% fp - fraction plume water
% zf - depths at which fs and fp values are calculated
% zs - depths of source shelf waters for mixing model
% zf - depths of source PMW waters for mixing model
% xf, xs and xp - indices of source fjord, shelf and plume waters for mixing model (i.e. inputs to
% mixing model are S_fjord(xf), S_shelf(xs), S_plume(xp), T_fjord(xf),
% T_shelf(xs), T_plume(xp))
% resn - sum of squared residuals from mixing model
% resi - individual residuals from mixing model

function[fs,fp,zf,zs,zp,xf,xs,xp,resn,resi] = mixing_model(S_fjord,T_fjord,rho0_fjord,z_fjord,S_shelf,T_shelf,rho0_shelf,z_shelf,S_plume,T_plume,rho0_plume,z_plume,z_fractions)
                                                            
%% find points of equal density

count = 1;

for z = z_fractions
    % xs and xp are vectors of cell locations in the shelf and plume profiles that have equal density to
    % the fjord profile, at selected intervals

    xf(count) = find(round(z_fjord)==round(z));
    zf(count) = z_fjord(xf(count)); % depth at which fjord water mass is found

    rhodiff = abs(rho0_shelf-rho0_fjord(xf(count)));
    xs(count) = find(rhodiff==min(rhodiff),1);
    zs(count) = z_shelf(xs(count)); % depth at which shelf water mass is found

    % only match with plume profile where lies within neutral buoyancy depth range and fjord density lies within bounds
    % of plume density
    if z > min(z_plume) && z < max(z_plume) && ...
            rho0_fjord(xf(count)) > min(rho0_plume) && rho0_fjord(xf(count)) < max(rho0_plume)
        rhodiff = abs(rho0_plume-rho0_fjord(xf(count)));
        xp(count) = find(rhodiff==min(rhodiff),1);
        zp(count) =  z_plume(xp(count)); % depth at which plume water mass is found
    else
        xp(count) = NaN;
        zp(count) = NaN;
    end

    count = count+1;
end

%% calculate fractions

for j = 1:length(xf)

    % get T values for water parcels of equal density
    T1 = T_fjord(xf(j)); % fjord
    T2 = T_shelf(xs(j)); % shelf
    if ~isnan(xp(j)) % if plume water present at this density
        T4 = T_plume(xp(j)); % plume
    else
        T4 = NaN; % plume
    end

    % get S values for water parcels of equal density
    S1 = S_fjord(xf(j)); % fjord
    S2 = S_shelf(xs(j)); % shelf
    if ~isnan(xp(j)) % if plume water present at this density
        S4 = S_plume(xp(j)); % plume
    else
        S4 = NaN; % plume
    end

    if ~isnan(T4) 
       
        A = [T2,T4;S2,S4;1,1];
        B = [T1;S1;1];
        
        [X,resn(j),resi(j,:)] = lsqlin(A,B,-ones(size(A)),zeros(size(B)),[],[],[0 0],[1 1]); % limit fm and fs to lie between 0 and 1, and fs+fm = 1

        fs(j) = X(1); % fraction shelf water
        fp(j) = X(2); % fraction plume water
        
    else % if no plume water at this density

        resn(j) = 0;
        resi(j,:) = [0 0 0];

        fs(j) = 1; % fraction shelf water
        fp(j) = 0; % fraction plume water
        
    end
end
