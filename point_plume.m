function diff=point_plume(z,X,alpha,g,rho0,Cd,GT,GS,cw,ci,lambda1,lambda2,lambda3,Tice,L,RHOa,Ta,Sa,depth)
% Donald Slater 02/14, updated Tom Cowton 03/23

diff=zeros(4,1);

% FIND Sb, Tb AND mdot AS DESCRIBED IN PDF
a=lambda1*(GT*cw-GS*ci);
b=GS*ci*(lambda1*X(4)-lambda2-lambda3*z+Tice-(L/ci))-GT*cw*(X(3)-lambda2-lambda3*z);
c=GS*X(4)*(ci*(lambda2+lambda3*z-Tice)+L);
Sb=(1/(2*a))*(-b-sqrt(b^2-4*a*c));
Tb=lambda1*Sb+lambda2+lambda3*z;
mdot=GS*sqrt(Cd)*X(2)*(X(4)-Sb)/Sb;

P = 1027*g*abs(z)*1*10^(-4); % pressure in dbar

% DEFINE MAIN EQUATIONS
diff(1)=2*alpha+4*mdot/(pi*X(2))-X(1)*g*(interp1(depth,RHOa,z)-gsw_rho(X(4),X(3),P))/(2*X(2)*X(2)*rho0)+2*Cd/pi;
diff(2)=-2*alpha*X(2)/X(1)-4*mdot/(pi*X(1))+g*(interp1(depth,RHOa,z)-gsw_rho(X(4),X(3),P))/(X(2)*rho0)-4*Cd*X(2)/(pi*X(1));
diff(3)=2*alpha*(interp1(depth,Ta,z)-X(3))/X(1)+4*mdot*(Tb-X(3))/(pi*X(1)*X(2))-4*GT*sqrt(Cd)*(X(3)-Tb)/(pi*X(1));
diff(4)=2*alpha*(interp1(depth,Sa,z)-X(4))/X(1)+4*mdot*(Sb-X(4))/(pi*X(1)*X(2))-4*GS*sqrt(Cd)*(X(4)-Sb)/(pi*X(1));