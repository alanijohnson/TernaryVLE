clear all;
clc;

%x and y values
x_acetone = linspace(0,1,100);
x_ethanol = linspace(0,1,100);
x_ethyl = linspace(0,1,100);
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
A12 = 0.9561; %acetone,ethanol
A21 = 1.6029; %ethanol,acetone
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
%system properties
P = 1.013; %bar
T = 298; %kelvin

%bubble point
for x_a = x_acetone
    for x_e = x_ethanol 
        for x_y = x_ethyl
            if (x_a+x_e+x_y==1)
             % temperature guess
             i=i+1;
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
                        T_liq(i) = T;
                        y_ace(i) = y_a;
                        y_eth(i) = y_e;
                        y_yl(i) = y_y;
                        liq_ace(i) = x_a;
                        liq_eth(i) = x_e;
                        liq_yl(i) = x_y;
                        T = T;
                        break;
                    end
                end
            end
        end
    end
end


%dew point
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
                    
%                     g1 = exp((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2);
%                     g2 = exp((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2);
%                     g3 = exp((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2);
%                     
                    x_ac = x_a*P/g1/P_vap_A; %x_a(x) = vapor mol fraction
                    x_et = x_e*P/g2/P_vap_E;
                    x_ey = (1-x_a-x_e)*P/g3/P_vap_Y;
                    sum_x = x_ac+x_et+x_ey;
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
                        T_vap(i) = T;
%                         x_ace(i) = x_ac;
%                         x_eth(i) = x_et;
%                         x_yl(i) = x_ey;
                          x_ace(i) = x_a;
                          x_eth(i) = x_e;
                          x_yl(i) = x_y;
                  
                        T = T;
                        break;
                    end
                end
            end
        end
    end
end



 disp(min(T_liq))
 disp(max(T_liq))
% b = true;
% val=110000
% while (b)
%     [val,ind]= max(T_liq(T_liq~=val));
%     if (x_ace(ind)~=0 && x_eth(ind)~=0 && (1-x_ace(ind)-x_eth(ind))~=0)
%         b=false;
%     end
% end

% disp(x_ace(ind))
% disp(x_eth(ind))
% disp(1-x_ace(ind)-x_eth(ind))
% disp(T_liq(ind))
% % disp(min(T_vap))
% disp(max(T_vap))

% %ternary plot Temp Contour Lines
% figure()
% terncontour(x_ace,x_eth,x_yl,T_liq)
% ternlabel('x_A_c_e_t_o_n_e', 'x_E_t_h_a_n_o_l', 'x_E_t_h_y_l_ _A_c_e_t_a_t_e')
% title('Temperature surface contours')
% %axis([0 1 0 1 0 1 0 1])

% %Bubble point Dew Point tern plot with color
% figure()
% ternpcolor(liq_ace, liq_eth, liq_yl,T_liq);
% shading interp
% ternlabel('x_A_c_e_t_o_n_e', 'x_E_t_h_a_n_o_l', 'x_E_t_h_y_l_ _A_c_e_t_a_t_e')
% title('Bubble Points for Acetone/Ethanol/Ethyl Acetate Mixtures')
% figure()
% ternpcolor(y_ace, y_eth, y_yl,T_vap);
% shading interp
% ternlabel('x_A_c_e_t_o_n_e', 'x_E_t_h_a_n_o_l', 'x_E_t_h_y_l_ _A_c_e_t_a_t_e')
% title('Dew Points for Acetone/Ethanol/Ethyl Acetate Mixtures')

figure()
ternsurf(liq_ace,liq_eth,liq_yl,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1 0 1])
legend('(T-329.418)/(351.4410-329.418)')
ternlabel2('x_A_c_e_t_o_n_e', 'x_E_t_h_a_n_o_l', 'x_E_t_h_y_l_ _A_c_e_t_a_t_e')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')

figure()
ternsurf(liq_eth,liq_ace,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1])
legend('(T-328.206)/(351.4410-328.206)')
ternlabel2('x_E_t_h_a_n_o_l','x_A_c_e_t_o_n_e',  'x_E_t_h_y_l_ _A_c_e_t_a_t_e')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')

figure()
ternsurf(1-liq_ace-liq_eth,liq_eth,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1])
legend('(T-328.206)/(351.4410-328.206)')
ternlabel2(  'x_E_t_h_y_l_ _A_c_e_t_a_t_e','x_E_t_h_a_n_o_l','x_A_c_e_t_o_n_e')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')

figure()
ternsurf(1-liq_ace-liq_eth,liq_ace,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1])
legend('(T-328.206)/(351.4410-328.206)')
ternlabel2('x_E_t_h_y_l_ _A_c_e_t_a_t_e','x_A_c_e_t_o_n_e','x_E_t_h_a_n_o_l')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')

figure()
ternsurf(liq_ace,1-liq_ace-liq_eth,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1])
legend('(T-328.206)/(351.4410-328.206)')
ternlabel2('x_A_c_e_t_o_n_e','x_E_t_h_y_l_ _A_c_e_t_a_t_e','x_E_t_h_a_n_o_l')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')

figure()
ternsurf(liq_eth,1-liq_ace-liq_eth,((T_liq-328.206)./(351.4410-328.206)))
shading interp
%axis([0 1 0 1])
legend('(T-328.206)/(351.4410-328.206)')
ternlabel2('x_E_t_h_a_n_o_l','x_E_t_h_y_l_ _A_c_e_t_a_t_e','x_A_c_e_t_o_n_e')
title('Temperature Surface Contour for Acetone/Ethanol/Ethyl Acetate Mixture')
% 

% 
% Y = array2table(y_a')
% Z = array2table(x_a')
% sprintf('%10.5f',T_liq)
%             
            
    
    
    
    
    