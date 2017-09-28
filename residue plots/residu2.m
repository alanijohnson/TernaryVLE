% Residue Curve Map for Ternary Systems
% Case of Acetone - Chloroform - Methanol System
% Author's Data: Housam BINOUS
% Department of Chemical Engineering
% National Institute of Applied Sciences and Technology
% Tunis, TUNISIA
% Email: binoushousam@yahoo.com 

%Acetone = component 1
%Chloroform = component 2
%Methanol = component 3 

function xdot=residu2(t,x,flag)


% Gas Constant and Total Pressure
R = 1.987; P = 760;

% Liquid mole fractions for Acetone - Chloroform - Methanol
% and Temperature are x(1), x(2), x(3) and x(4)
T=x(4);

if isempty(flag),
    
% Vapor Pressure Using Antoine Equation    

A1 = 7.11714; B1 = 1210.595; C1 = 229.664;
A2 = 6.95465; B2 = 1170.966; C2 = 226.232;
A3 = 8.08097; B3 = 1582.271; C3 = 239.726;

Psat1=10.^(A1-B1/(C1+T));
Psat2=10.^(A2-B2/(C2+T));
Psat3=10.^(A3-B3/(C3+T));

% Binary inetraction Parameters and Activities for Wilson Model

l12 = 116.1171; l21 = -506.8519; l13 = -114.4047;
l31 = 545.2942; l23 = -361.7944; l32 = 1694.0241;

V1 = 74.05; V2 = 80.67; V3 = 40.73;

A12 = V2/V1*exp(-l12/(R*(273.15 + T)));
A13 = V3/V1*exp(-l13/(R*(273.15 + T)));
A32 = V2/V3*exp(-l32/(R*(273.15 + T)));
A21 = V1/V2*exp(-l21/(R*(273.15 + T)));
A31 = V1/V3*exp(-l31/(R*(273.15 + T)));
A23 = V3/V2*exp(-l23/(R*(273.15 + T)));

gamma1 = exp(-log(x(1) + x(2)*A12 + x(3)*A13) + 1 - (x(1)/(x(1) + x(2)*A12 + x(3)*A13) + ... 
x(2)*A21/(x(1)*A21 + x(2) + x(3)*A23) + x(3)*A31/(x(1)*A31 + x(2)*A32 + x(3))));

gamma2 = exp(-log(x(1)*A21 + x(2) + x(3)*A23) + 1 - (x(1)*A12/(x(1) + x(2)*A12 + x(3)*A13) + ...
x(2)/(x(1)*A21 + x(2) + x(3)*A23) + x(3)*A32/(x(1)*A31 + x(2)*A32 + x(3))));

gamma3 = exp(-log(x(1)*A31 + x(2)*A32 + x(3)) + 1 - (x(1)*A13/(x(1) + x(2)*A12 + x(3)*A13) + ...
x(2)*A23/(x(1)*A21 + x(2) + x(3)*A23) + x(3)/(x(1)*A31 + x(2)*A32 + x(3))));

% Vapor mole fractions for Acetone - Chloroform - Methanol

y(1)=x(1)*gamma1*Psat1/P;
y(2)=x(2)*gamma2*Psat2/P;
y(3)=x(3)*gamma3*Psat3/P;

% Differential Algeraic Equations for residue curve map calculation

xdot(1)=y(1)-x(1);
xdot(2)=y(2)-x(2);
xdot(3)=y(3)-x(3);
xdot(4) = x(1)*gamma1*Psat1 + x(2)*gamma2*Psat2 + x(3)*gamma3*Psat3 - P;

xdot = xdot';  % xdot must be a column vector

else
% Return M

   M = zeros(4,4);
   for i = 1:3,
      M(i,i) = 1;
   end
   xdot = M;
end
