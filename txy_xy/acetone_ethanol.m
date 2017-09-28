clear all;
clc;

%x and y values
x_a = linspace(0,.05,20);
T_vap = [298];
T_liq = [298];


%acetone properties - component 1
A_A = 4.42448;
A_B = 1312.253;
A_C = -32.445;

%ethanol properties - component 2
E_A = 5.37229;
E_B = 1670.409;
E_C = -40.191;

%Van laar parameters
A12 = 0.5615;
A21 = 0.6580;

%Gamma values
G1 = exp(A12.*(A21.*(1-x_a)./(A12.*x_a+A21.*(1-x_a))).^2);
G2 = exp(A21.*(A12.*(x_a)./(A12.*x_a+A21.*(1-x_a))).^2);

%system properties
P = 1; %1 atm
T = 298;
%index the x vals
for x=1:length(x_a)
     % temperature guess
     i = length(T)-1;
    while(1)
        i = i+1;
        %disp(i);
        P_vap_A = 10^(A_A-(A_B/(A_C+T)));
        P_vap_E = 10^(E_A-(E_B/(E_C+T)));
        g1 = exp(A12*(A21*(1-x_a(x))/(A12*x_a(x)+A21*(1-x_a(x))))^2);
        g2 = exp(A21*(A12*(x_a(x))/(A12*x_a(x)+A21*(1-x_a(x))))^2);
        y_ace = x_a(x)*g1*P_vap_A/P;
        y_e = (1-x_a(x))*g2*P_vap_E/P;
        sum_y = y_ace+y_e;
        if sum_y>1.1
            T = T-0.1;
        elseif sum_y<0.9
            T = T+0.1;
        elseif sum_y>1.01
            T = T-0.01;
        elseif sum_y<0.99
            T = T + 0.01;
        elseif sum_y<0.99999
            T = T + 0.00001;
        elseif sum_y>1.00001
            T = T-0.00001;
        else
            T_vap(x) = T;
            y_a(x) = y_ace;
            T = T;
            break;
        end
    end
end


T = 298;
i = length(T)-1;
%index the x vals
for x=1:length(x_a)
    i = length(T)-1;
    while(1)
        i=i+1;
        %disp(i)
        P_vap_A = 10^(A_A-(A_B/(A_C+T)));
        P_vap_E = 10^(E_A-(E_B/(E_C+T)));
        g1 = exp(A12*(A21*(1-x_a(x))/(A12*x_a(x)+A21*(1-x_a(x))))^2);
        g2 = exp(A21*(A12*(x_a(x))/(A12*x_a(x)+A21*(1-x_a(x))))^2);
        x_ace = x_a(x)*P/g1/P_vap_A; %x_a(x) = vapor mol fraction
        x_e = (1-x_a(x))*P/g2/P_vap_E;
        sum_x = x_ace+x_e;
        if sum_x> 1.1
           T = T+0.1;
        elseif sum_x< 0.9
           T = T-0.1;
        elseif sum_x>1.01
            T = T+0.01;
        elseif sum_x<0.99
            T = T-0.01;
        elseif sum_x>1.00001
            T = T+0.00001;
        elseif sum_x<0.99999
            T = T-0.00001;
        else
            T_liq(x) = T;
            T = T+50;
            break;
        end
    end
end

figure()
plot(x_a,T_vap,'--')
hold on
plot(x_a,T_liq)

title('T-x-y for Acetone(1)-Ethanol(2) Mixture')
xlabel('x1/y1')
ylabel('Temperature (K)')
legend('liquid','vapor')

figure()
plot(x_a,y_a)
hold on
plot (x_a,x_a,'--')
xlabel('x1')
ylabel('y1')
title('x-y for Acetone(1)-Ethanol(2) Mixture')
legend('x/y','x=y','Location','southeast')
axis([0 1 0 1])

figure()
plot(x_a,G1)
hold on
plot(x_a,G2)
xlabel('x1')
ylabel('gamma')
title('Activity Coefficients for Acetone(1)-Ethanol(2) Mixture')
legend('gamma 1','gamma 2','Location','northeast')
axis([0 1 1 2])

Y = array2table(x_a')
Z = array2table(T_liq')
sprintf('%10.5f',T_liq)
            
            
    
    
    
    
    