function diff=line_plume(z,X,RHO,T,S,depth,E_0,k,rho_ref,g,cw,ci,L,GT,GS,lambda1,lambda2,lambda3,Tice)
% Donald Slater 02/14, updated Tom Cowton 03/23

diff=zeros(4,1);

% FIND Sb, Tb AND mdot AS DESCRIBED IN PDF
a=lambda1*(GT*cw-GS*ci);
b=GS*ci*(lambda1*X(4)-lambda2-lambda3*z+Tice-(L/ci))-GT*cw*(X(3)-lambda2-lambda3*z);
c=GS*X(4)*(ci*(lambda2+lambda3*z-Tice)+L);
Sb=(1/(2*a))*(-b-sqrt(b^2-4*a*c));
Tb=lambda1*Sb+lambda2+lambda3*z;
mdot=GS*sqrt(k)*X(2)*(X(4)-Sb)/Sb;

P = 1027*g*abs(z)*1*10^(-4); % pressure in dbar

% DEFINE EQUATIONS
diff(1)=2*E_0+k-(g*X(1)/(X(2)^2))*(interp1(depth,RHO,z)-gsw_rho(X(4),X(3),P))/rho_ref+2*mdot/X(2);
diff(2)=-(X(2)/X(1))*(E_0+k+mdot/X(2))+(g/X(2))*(interp1(depth,RHO,z)-gsw_rho(X(4),X(3),P))/rho_ref;
diff(3)=E_0*interp1(depth,T,z)/X(1)-(X(3)/X(1))*(E_0+mdot/X(2))+(mdot/(X(1)*X(2)))*(Tb-(L/cw)-(ci/cw)*(Tb-Tice));
diff(4)=E_0*interp1(depth,S,z)/X(1)-(X(4)/X(1))*(E_0+mdot/X(2));