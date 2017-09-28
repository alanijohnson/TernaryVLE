x_e = 0.9
x_a = 0.01
x_y = 0.09

P=1.013
T=350
while(1)
                    %disp(i)
                    P_vap_A = 10^(A_A-(A_B/(A_C+T)));
                    P_vap_E = 10^(E_A-(E_B/(E_C+T)));
                    P_vap_Y = 10^(Y_A-(Y_B/(Y_C+T)));
                    %((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+((x_e*x_y*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+x_y*(A31/A13))^2)
                    
                    g1 = exp(((x_e^2*A12*(A21/A12)^2)+(x_y^2*A13*(A31/A13)^2)+(x_e*x_y*A21*A31/A12/A13*(A12+A13-A23*A12/A21)))/((x_a+x_e*A21/A12+x_y*A31/A13)^2))
                    g2 = exp(((x_a^2*A21*(A12/A21)^2)+(x_y^2*A23*(A32/A23)^2)+(x_a*x_y*A12*A32/A21/A23*(A21+A23-A13*A21/A12)))/((x_e+x_a*A12/A21+x_y*A32/A23)^2))
                    g3 = exp(((x_y^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+(x_e*x_a*A23*A13/A32/A31*(A32+A23-A21*A32/A23)))/((x_y+x_e*A23/A32+x_a*A13/A31)^2))
                    
%                     g1 = exp((x_e^2*A12*(A21/A12)^2)+((1-x_a-x_e)^2*A13*(A31/A13)^2)+((x_e*(1-x_a-x_e)*A21*A31/A12/A13)*(A12+A13-A23*A12/A21)))/((x_a+x_e*(A21/A12)+(1-x_a-x_e)*(A31/A13))^2);
%                     g2 = exp((x_a^2*A21*(A12/A21)^2)+((1-x_a-x_e)^2*A23*(A32/A23)^2)+((x_a*(1-x_a-x_e)*A12*A32/A21/A23)*(A21+A23-A13*A21/A12)))/((x_e+x_a*(A12/A21)+(1-x_a-x_e)*(A32/A23))^2);
%                     g3 = exp((x_e^2*A32*(A23/A32)^2)+(x_a^2*A31*(A13/A31)^2)+((x_e*x_a*A23*A13/A32/A31)*(A32+A31-A21*A32/A23)))/(((1-x_a-x_e)+x_e*(A23/A32)+x_a*(A13/A31))^2);
%                     
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
                        T_liq(i) = T;
                        x_ace(i) = x_ac;
                        x_eth(i) = x_et
                        T = T;
                        break;
                    end
                end