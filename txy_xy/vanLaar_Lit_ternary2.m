clear all;
clc;

%x and y values
x_acetone = linspace(0,1,1000);
x_ethanol = linspace(0,1,1000);
x_ethyl = linspace(0,1,1000);
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

%ethyl acetate properties - component 3
Y_A = 4.22809;
Y_B = 1245.702;
Y_C = -55.189;

%Van laar parameters
A12 = 0.5615; %acetone,ethanol
A21 = 0.6580; %ethanol,acetone
A13 = 0.1919; %acetone,ethyl acetate
A31 = 0.161; %ethyl acetate, acetone
A23 = 0.565; %ethanol, ethyl acetate
A32 = 0.913; %ethyl acetate, ethanol

% %Gamma values
%         
% 
% G1 = exp(A12.*(A21.*(1-x_a)./(A12.*x_a+A21.*(1-x_a))).^2);
% G2 = exp(A21.*(A12.*(x_a)./(A12.*x_a+A21.*(1-x_a))).^2);
% G3 =
i=0;
j=0;
k=0;
l=0;
m=0;
n=0;
o=0;
p=0;
q=0;
A_y_ey=0;
A_y_ey95=0;
A_y_ey105=0;
B_y_ey=0;
B_y_ey95=0;
B_y_ey105=0;
C_y_ey=0;
C_y_ey95=0;
C_y_ey105=0;
A_y_ace=0;
A_y_ace95=0;
A_y_ace105=0;
A_y_eth=0;
A_y_eth95=0;
A_y_eth105=0;

B_y_ace=0;
B_y_ace95=0;
B_y_ace105=0;
B_y_eth=0;
B_y_eth95=0;
B_y_eth105=0;

C_y_ace=0;
C_y_ace95=0;
C_y_ace105=0;
C_y_eth=0;
C_y_eth95=0;
C_y_eth105=0;
%system properties
P = 1.013; %bar
T = 298; %kelvin

%dew point
for x_a = x_acetone
    for x_e = x_ethanol
        for x_y = x_ethyl
            if (x_a+x_e+x_y==1)
                % temperature guess
                
                while(1)
                    %i = i+1;
                    %disp(i);
                    P_vap_A = 10^(A_A-(A_B/(A_C+T)));
                    P_vap_E = 10^(E_A-(E_B/(E_C+T)));
                    P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
                    %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
                    g1 = exp(((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2));
                    g2 = exp(((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2));
                    g3 = exp(((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2));
                    
                    y_a = x_a*g1*P_vap_A/P;
                    y_e = (x_e)*g2*P_vap_E/P;
                    y_y = (1-x_a-x_e)*g3*P_vap_Y/P;
                    sum_y = y_a+y_e+y_y;
                    if sum_y>1.1
                        T = T-1;
                    elseif sum_y<0.9
                        T = T+1;
                    elseif sum_y>1.01
                        T = T-0.1;
                    elseif sum_y<0.99
                        T = T + 0.1;
                    elseif sum_y<0.999
                        T = T + 0.01;
                    elseif sum_y>1.001
                        T = T-0.01;
                    elseif sum_y<0.99975
                        T = T + 0.001;
                    elseif sum_y>1.00025
                        T = T-0.001;
                    else
                        k1 = (g1*P_vap_A/P);
                        k2 = (g2*P_vap_E/P);
                        k3 = (g3*P_vap_Y/P);
                        %k=1
                        if (abs(1-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            i=i+1;
                            A_y_ace(i)=x_a;
                            A_y_eth(i)=x_e;
                            A_y_ey(i)= x_y;
                        elseif (abs(1-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            j=j+1;
                            B_y_ace(j)=x_a;
                            B_y_eth(j)=x_e;
                            B_y_ey(j)=x_y;
                        elseif (abs(1-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            k=k+1;
                            C_y_ace(k)=x_a;
                            C_y_eth(k)=x_e;
                            C_y_ey(k) = x_y;
                        %0.95 limits
                        elseif (abs(0.95-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            l=l+1;
                            A_y_ace95(l)=x_a;
                            A_y_eth95(l)=x_e;
                            A_y_ey95(l)=x_y;
                        elseif (abs(0.95-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            m=m+1;
                            B_y_ace95(m)=x_a;
                            B_y_eth95(m)=x_e;
                            B_y_ey95(m)=x_y;
                        elseif (abs(0.95-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            n=n+1;
                            C_y_ace95(n)=x_a;
                            C_y_eth95(n)=x_e;
                            C_y_ey95(n)=x_y;
                        %1.05 limits
                        elseif (abs(1.05-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            o=o+1;
                            A_y_ace105(o)=x_a;
                            A_y_eth105(o)=x_e;
                            A_y_ey105(o)=x_y
                        elseif (abs(1.05-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            p=p+1;
                            B_y_ace105(p)=x_a;
                            B_y_eth105(p)=x_e;
                            B_y_ey105(p)=x_y;
                        elseif (abs(1.05-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            q=q+1;
                            C_y_ace105(q)=x_a;
                            C_y_eth105(q)=x_e;
                            C_y_ey105(q)=x_y;
                        end
                        T = T;
                        break;
                    end
                end
            end
        end
    end
end

%bubble point
T = 298;
i = 0;
%index the x vals
for x_a = x_acetone
    for x_e = x_ethanol
        for x_y= x_ethyl
            if (x_a+x_e+x_y==1)
                i = i+1;
                while(1)
                    %disp(i)
                    P_vap_A = 10^(A_A-(A_B/(A_C+T)));
                    P_vap_E = 10^(E_A-(E_B/(E_C+T)));
                    P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
                    %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
                    g1 = exp(((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2));
                    g2 = exp(((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2));
                    g3 = exp(((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2));
                    
                    x_ac = x_a*P/g1/P_vap_A; %x_a(x) = vapor mol fraction
                    x_et = x_e*P/g2/P_vap_E;
                    x_y = (1-x_a-x_e)*P/g3/P_vap_Y;
                    sum_x = x_ac+x_et+x_y;
                    if sum_x> 1.1
                        T = T+1;
                    elseif sum_x< 0.9
                        T = T-1;
                    elseif sum_x>1.01
                        T = T+0.1;
                    elseif sum_x<0.99
                        T = T-0.1;
                    elseif sum_x>1.001
                        T = T+0.01;
                    elseif sum_x<0.999
                        T = T-0.01;
                    elseif sum_x>1.00025
                        T=T+0.001;
                    elseif sum_x<0.99975
                        T=T-0.001;
                    else
                        k1 = (g1*P_vap_A/P);
                        k2 = (g2*P_vap_E/P);
                        k3 = (g3*P_vap_Y/P);
                        %k=1
                        if (abs(1-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            i=i+1;
                            A_y_ace(i)=x_a;
                            A_y_eth(i)=x_e;
                            A_y_ey(i)= x_y;
                        elseif (abs(1-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            j=j+1;
                            B_y_ace(j)=x_a;
                            B_y_eth(j)=x_e;
                            B_y_ey(j)=x_y;
                        elseif (abs(1-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            k=k+1;
                            C_y_ace(k)=x_a;
                            C_y_eth(k)=x_e;
                            C_y_ey(k) = x_y;
                        %0.95 limits
                        elseif (abs(0.95-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            l=l+1;
                            A_y_ace95(l)=x_a;
                            A_y_eth95(l)=x_e;
                            A_y_ey95(l)=x_y;
                        elseif (abs(0.95-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            m=m+1;
                            B_y_ace95(m)=x_a;
                            B_y_eth95(m)=x_e;
                            B_y_ey95(m)=x_y;
                        elseif (abs(0.95-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            n=n+1;
                            C_y_ace95(n)=x_a;
                            C_y_eth95(n)=x_e;
                            C_y_ey95(n)=x_y;
                        %1.05 limits
                        elseif (abs(1.05-(g1*P_vap_A/P))<0.001) %k=1 for acetone
                            o=o+1;
                            A_y_ace105(o)=x_a;
                            A_y_eth105(o)=x_e;
                            A_y_ey105(o)=x_y
                        elseif (abs(1.05-(g2*P_vap_E/P))<0.001) %k=1 for ethanol
                            p=p+1;
                            B_y_ace105(p)=x_a;
                            B_y_eth105(p)=x_e;
                            B_y_ey105(p)=x_y;
                        elseif (abs(1.05-(g3*P_vap_Y/P))<0.001) %k=1 for ethyl acetate
                            q=q+1;
                            C_y_ace105(q)=x_a;
                            C_y_eth105(q)=x_e;
                            C_y_ey105(q)=x_y;
                        end
                        T = T;
                        break;
                    end
                end
            end
        end
    end
end




%Acetone K=1

figure()
ternplot(A_y_ace,A_y_eth,A_y_ey,'g.')
hold on
ternplot(A_y_ace95,A_y_eth95,A_y_ey95,'r.')
hold on
ternplot(A_y_ace105,A_y_eth95,A_y_ey105,'r.')
ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
title('Acetone K=1 in Acetone/Ethanol/Ethyl Acetate Mixture')
legend('K=1','K=0.95','K=1.05')


%Ethanol K=1


figure()
ternplot(B_y_ace,B_y_eth,B_y_ey,'g.')
hold on
ternplot(B_y_ace95,B_y_eth95,B_y_ey95,'r.')
hold on
ternplot(B_y_ace105,B_y_eth105,B_y_ey105,'r.')
ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
title('Ethanol K=1 in Acetone/Ethanol/Ethyl Acetate Mixture')
legend('K=1','K=0.95','K=1.05')

%Ethyl Acetate K=1
figure()
ternplot(C_y_ace,C_y_eth,C_y_ey,'g.')
hold on
ternplot(C_y_ace95,C_y_eth95,C_y_ey95,'r.')
hold on
ternplot(C_y_ace105,C_y_eth105,C_y_ey105,'r.')
ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
title('Ethyl Acetate K=1 in Acetone/Ethanol/Ethyl Acetate Mixture')
legend('K=1','K=0.95','K=1.05')


%All on one plot
figure()
ternplot(A_y_ace,A_y_eth,1-A_y_ace-A_y_eth,'r.')
hold on
ternplot(B_y_ace,B_y_eth,1-B_y_ace-B_y_eth,'g.')
hold on
ternplot(C_y_ace,C_y_eth,1-C_y_ace-C_y_eth,'b.')
ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
legend('K=1 Acetone', 'K=1 Ethanol', 'K=1 Ethyl Acetate')
title('K=1 for components in Acetone/Ethanol/Ethyl Acetate Mixture')
legend('K=1','K=0.95','K=1.05')





% %bubble point
% T = 298;
% i = 0;
% %index the x vals
% for x_a = x_acetone
%     for x_e = x_ethanol
%         for x_y= x_ethyl
%             if (x_a+x_e+x_y==1)
%                 i = i+1;
%                 while(1)
%                     %disp(i)
%                     P_vap_A = 10^(A_A-(A_B/(A_C+T)));
%                     P_vap_E = 10^(E_A-(E_B/(E_C+T)));
%                     P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
%                     %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
%                     g1 = exp((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2);
%                     g2 = exp((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2);
%                     g3 = exp((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2);
%                     
%                     x_ac = x_a*P/g1/P_vap_A; %x_a(x) = vapor mol fraction
%                     x_et = x_e*P/g2/P_vap_E;
%                     x_y = (1-x_a-x_e)*P/g3/P_vap_Y;
%                     sum_x = x_ac+x_et+x_y;
%                     if sum_x> 1.1
%                         T = T+1;
%                     elseif sum_x< 0.9
%                         T = T-1;
%                     elseif sum_x>1.01
%                         T = T+0.1;
%                     elseif sum_x<0.99
%                         T = T-0.1;
%                     elseif sum_x>1.001
%                         T = T+0.01;
%                     elseif sum_x<0.999
%                         T = T-0.01;
%                     elseif sum_x>1.00025
%                         T=T+0.001;
%                     elseif sum_x<0.99975
%                         T=T-0.001;
%                     else
%                         T_liq(i) = T;
%                         x_ace(i) = x_ac;
%                         x_eth(i) = x_et;
%                   
%                         T = T;
%                         break;
%                     end
%                 end
%             end
%         end
%     end
% end



% %contour3(x_ace,x_eth,T_liq)
% figure()
% terncontour(x_ace,x_eth,(1-x_ace-x_eth),T_liq)
% ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
% title('Temperature surface contours')
% axis([0 1])
% figure()
% ternpcolor(x_ace, x_eth, (1-x_ace-x_eth),T_liq);
% shading interp
% ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
% title('Bubble Points for Acetone/Ethanol/Ethyl Acetate Mixtures')
% figure()
% ternpcolor(y_ace, y_eth, (1-x_ace-x_eth),T_vap);
% shading interp
% ternlabel('Acetone', 'Ethanol', 'Ethyl Acetate')
% title('Dew Points for Acetone/Ethanol/Ethyl Acetate Mixtures')
% figure()
% ternsurf(x_ace,x_eth,T_liq)
% axis([0 1 0 1])
% % figure()
% % plot(x,T_vap,'--')
% % hold on
% plot(x_a,T_liq)
% 
% title('T-x-y for Acetone(1)-Ethanol(2) Mixture')
% xlabel('x1/y1')
% ylabel('Temperature (K)')
% legend('liquid','vapor')
% 
% figure()
% plot(x_a,y_ace)
% hold on
% plot (x_a,x_a,'--')
% xlabel('x1')
% ylabel('y1')
% title('x-y for Acetone(1)-Ethanol(2) Mixture')
% legend('x/y','x=y','Location','southeast')
% axis([0 1 0 1])

% figure()
% plot(x_a,G1)
% hold on
% plot(x_a,G2)
% xlabel('x1')
% ylabel('gamma')
% title('Activity Coefficients for Acetone(1)-Ethanol(2) Mixture')
% legend('gamma 1','gamma 2','Location','northeast')
% axis([0 1 1 2])

Y = array2table(x_a')
Z = array2table(T_liq')
sprintf('%10.5f',T_liq)
            
            
    
    
    
    
    