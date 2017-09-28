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

% Main program calls residu and residu2 and uses ode15s to solve the DAE.
% Different initial conditions are taken in order to plot several residu curves

for i=0:9,
xO1=i*0.1;
x02=0.1*(9-i);
x03=1-0.1*i-0.1*(9-i);
T0=57.3123;

tf=10;

x0 = [xO1  x02  x03  T0];

opts = odeset('Mass','M','MassSingular','yes');

[t,x] = ode15s('residu',[0 tf],x0,opts);


x1=x(:,1);
x2=x(:,2);

figure(1);
x1=x(:,1);
x2=x(:,2);

AXIS([0 1 0 1])
hold on

plot(x1,x2,'r')
plot(line([0 1],[1 0]),'b');

tf=20

[t,x] = ode15s('residu2',[0 tf],x0,opts);
x1=x(:,1);
x2=x(:,2);
plot(x1,x2,'r')

end

for i=0:9,
    
xi3=i*0.1;
xi2=0.1*(9-i);
xi1=1-0.1*i-0.1*(9-i);
Ti=70;

tf=10;

x0 = [xi1  xi2  xi3  Ti];

opts = odeset('Mass','M','MassSingular','yes');

[t,x] = ode15s('residu',[0 tf],x0,opts);


x1=x(:,1);
x2=x(:,2);

figure(1);
x1=x(:,1);
x2=x(:,2);
hold on

plot(x1,x2,'r')
plot(line([0 1],[1 0]),'b');

tf=20

[t,x] = ode15s('residu2',[0 tf],x0,opts);
x1=x(:,1);
x2=x(:,2);
plot(x1,x2,'r')

end
for i=0:4,
    
xa2=(5-i)*0.1;
xa1=0.1*(2.5);
xa3=1-0.1*(5-i)-0.1*(2.5);
Ta=70;

tf=12.5;

x0 = [xa1  xa2  xa3  Ta];

opts = odeset('Mass','M','MassSingular','yes');

[t,x] = ode15s('residu',[0 tf],x0,opts);


x1=x(:,1);
x2=x(:,2);

figure(1);
x1=x(:,1);
x2=x(:,2);
hold on

plot(x1,x2,'r')
plot(line([0 1],[1 0]),'b');

tf=60

[t,x] = ode15s('residu2',[0 tf],x0,opts);
x1=x(:,1);
x2=x(:,2);
plot(x1,x2,'r')

end

hold off

