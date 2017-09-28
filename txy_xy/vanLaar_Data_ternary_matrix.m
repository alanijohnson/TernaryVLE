clear all;
clc;

%x and y values
x_a = linspace(0,1,10);
x_e = linspace(0,1,10);
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
A13 = 0.3611; %acetone,ethyl acetate
A31 = 0.0332; %ethyl acetate, acetone
A23 = 0.565; %ethanol, ethyl acetate
A32 = 0.913; %ethyl acetate, ethanol

% %Gamma values
%         
% 
g1 = exp((x_e.^2.*A12.*(A21/A12).^2)+((1-x_a-x_e).^2.*A13.*(A31/A13).^2)+((x_e.*(1-x_a-x_e).*A21.*A31./A12./A13).*(A12+A13-A23*A12/A21)))./((x_a+x_e.*(A21/A12)+(1-x_a-x_e).*(A31/A13)).^2);
g2 = exp((x_a.^2.*A21.*(A12/A21).^2)+((1-x_a-x_e).^2.*A23.*(A32/A23).^2)+((x_a.*(1-x_a-x_e).*A12.*A32./A21./A23).*(A21+A23-A13*A21/A12)))./((x_e+x_a.*(A12/A21)+(1-x_a-x_e).*(A32/A23)).^2);
g3 = exp((x_e.^2.*A32.*(A23/A32).^2)+(x_a.^2.*A31.*(A13/A31)^2)+((x_e.*x_a.*A23.*A13./A32./A31).*(A32+A31-A21*A32/A23)))./(((1-x_a-x_e)+x_e.*(A23/A32)+x_a.*(A13/A31)).^2);

surf(x_a,x_e,g1)
i=0;
j=0;
%system properties
P = 1.013; %bar
T = 298; %kelvin

% %dew point
% for x_a = x_acetone
%     i=i+1;
%     for x_e = x_ethanol 
%         if (x_a+x_e<=1)
%          % temperature guess
%          j=j+1;
%             while(1)
%                 %i = i+1;
%                 %disp(i);
%                 P_vap_A = 10^(A_A-(A_B/(A_C+T)));
%                 P_vap_E = 10^(E_A-(E_B/(E_C+T)));
%                 P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
%                 %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
%                 g1 = exp((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2);
%                 g2 = exp((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2);
%                 g3 = exp((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2);
% 
%                 y_a = x_a*g1*P_vap_A/P;
%                 y_e = (x_e)*g2*P_vap_E/P;
%                 y_y = (1-x_a-x_e)*g3*P_vap_Y/P;
%                 sum_y = y_a+y_e+y_y;
%                 if sum_y>1.1
%                     T = T-1;
%                 elseif sum_y<0.9
%                     T = T+1;
%                 elseif sum_y>1.01
%                     T = T-0.1;
%                 elseif sum_y<0.99
%                     T = T + 0.1;
%                 elseif sum_y<0.999
%                     T = T + 0.01;
%                 elseif sum_y>1.001
%                     T = T-0.01;
%                 elseif sum_y<0.99975
%                     T = T + 0.001;
%                 elseif sum_y>1.00025
%                     T = T-0.001;    
%                 else
%                     T_vap(i,j) = T;
%                     y_ace(i) = y_a;
%                     y_eth(j) = y_e;
%                     T = T;
%                     break;
%                 end
%             end
%         end
%     end
% end
% 
% 
% %bubble point
% T = 298;
% i = 0;
% j=0;
% %index the x vals
% for x_a = x_acetone
%     i=i+1;
%     for x_e = x_ethanol
%         j = j+1;
%         if (x_a+x_e<=1)
%             while(1)
%                 %disp(i)
%                 P_vap_A = 10^(A_A-(A_B/(A_C+T)));
%                 P_vap_E = 10^(E_A-(E_B/(E_C+T)));
%                 P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
%                 %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
%                 g1 = exp((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2);
%                 g2 = exp((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2);
%                 g3 = exp((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2);
% 
%                 x_ac = x_a*P/g1/P_vap_A; %x_a(x) = vapor mol fraction
%                 x_et = x_e*P/g2/P_vap_E;
%                 x_y = (1-x_a-x_e)*P/g3/P_vap_Y;
%                 sum_x = x_ac+x_et+x_y;
%                 if sum_x> 1.1
%                    T = T+1;
%                 elseif sum_x< 0.9
%                    T = T-1;
%                 elseif sum_x>1.01
%                     T = T+0.1;
%                 elseif sum_x<0.99
%                     T = T-0.1;
%                 elseif sum_x>1.001
%                     T = T+0.01;
%                 elseif sum_x<0.999
%                     T = T-0.01;
%                 elseif sum_x>1.00025
%                     T=T+0.001;
%                 elseif sum_x<0.99975
%                     T=T-0.001;
%                 else
%                     T_liq(i,j) = T;
%                     x_ace(i) = x_ac;
%                     x_eth(j) = x_et;
%                     T = T;
%                     break;
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% contour3(x_ace,x_eth,T_liq)
% %terncontour(x_ace,x_eth,(1-x_ace-x_eth),T_liq)
% % figure()
% % plot(x,T_vap,'--')
% % hold on
% % plot(x_a,T_liq)
% % 
% % title('T-x-y for Acetone(1)-Ethanol(2) Mixture')
% % xlabel('x1/y1')
% % ylabel('Temperature (K)')
% % legend('liquid','vapor')
% % 
% % figure()
% % plot(x_a,y_ace)
% % hold on
% % plot (x_a,x_a,'--')
% % xlabel('x1')
% % ylabel('y1')
% % title('x-y for Acetone(1)-Ethanol(2) Mixture')
% % legend('x/y','x=y','Location','southeast')
% % axis([0 1 0 1])
% 
% % figure()
% % plot(x_a,G1)
% % hold on
% % plot(x_a,G2)
% % xlabel('x1')
% % ylabel('gamma')
% % title('Activity Coefficients for Acetone(1)-Ethanol(2) Mixture')
% % legend('gamma 1','gamma 2','Location','northeast')
% % axis([0 1 1 2])
% 
% Y = array2table(x_a')
% Z = array2table(T_liq')
% sprintf('%10.5f',T_liq)
            
            
    
    
    
    
    