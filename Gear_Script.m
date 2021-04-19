%Strength and Stress Analysis for Set of Gears
%Author - Brandon Baeder

%Be sure to change gear properties (all found in lines 1-51)
%Be sure to change strength coefficients based on gear selection
%Strength coefficients are a,b,c,d (Lines 67-68, and 79-80)
%Final power rating is in horsepower
%Final surface and bending stresses are in psi



clear all;                                  %#ok<CLALL>
close all;
clc;

omega_nc = 1725;                            %rev/min
Power_nc = 0.5;                            %HP

S_ut_1 = 82000;                               %psi
S_e_in_1 = 0.5*S_ut_1;                            %psi (endurance limit set 1)   

S_ut_2 = 78000;                                 %psi
S_e_in_2 = 0.5*S_ut_2;                          %psi    (endurance limit set 2)

BH_1 = 122;                                   %Brinell Hardness Rating gearset 1
BH_2 = 130;                                   %brinell hardness rating gearset 2

omega = omega_nc*(2*pi);                    %rad/min
Power = 0.5*33000;                           %converted power to ft-lbf/min

phi = 20;                                   %degrees
n_s = 2;                                   %safety factor

N_1 = 12;                                   %number of teeth gear 1
F_1 = 0.125;                                %Face Width gear 1(in)
d_p_1 = 0.25;                               %Pitch Diameter gear 1 (in)
p_d_1 = N_1/d_p_1;                          %Diametral Pitch gear 1(Teeth/in)

N_2 = 48;                                   %number of teeth gear 2
F_2 = 0.25;                                 %Face Width gear 2(in)
d_p_2 = 0.92;                               %Pitch Diameter gear 2 (in)
p_d_2 = N_2/d_p_2;                          %Diametral Pitch gear 2(Teeth/in)

N_3 = 12;                                   %number of teeth gear 3
F_3 = 0.250;                                %Face Width gear 3(in)
d_p_3 = 0.50;                               %Pitch Diameter gear 3 (in)
p_d_3 = N_3/d_p_3;                          %Diametral Pitch gear 3(Teeth/in)

N_4 = 60;                                   %number of teeth gear 4
F_4 = 0.500;                                %Face Width gear 4(in)
d_p_4 = 3.00;                               %Pitch Diameter gear 4 (in)
p_d_4 = N_4/d_p_4;                         %Diametral Pitch gear 4(Teeth/in)

Y_1 = 0.245;                                %various lewis form factors
Y_2 = 0.403;
Y_3 = 0.245;
Y_4 = 0.422;

v_1 = omega*(d_p_1/24);                      %various linear velocities for gears (ft/min)           
v_2 = omega*(d_p_2/24);
v_3 = omega*(d_p_3/24);
v_4 = omega*(d_p_4/24);

Kv_1 = (50 + sqrt(v_1))/50;                 %various velocity factors 
Kv_2 = (50 + sqrt(v_2))/50;  
Kv_3 = (50 + sqrt(v_3))/50;  
Kv_4 = (50 + sqrt(v_4))/50;

W_t_1 = Power/v_1
W_t_3 = Power/v_2                     %various transmitted loads
W_t_2 = Power/v_3
W_t_4 = Power/v_4
a = 1.34;                                   %constants for strength analysis
b = -0.085;

k_a_1 = a*(S_ut_1^b);                           %fatigue coefficients       
k_b_1 = 1;
k_c_1 = 1;
k_d_1 = 1;
k_e_1 = 1;
k_f_1 = 1.89;

S_e_1= k_a_1*k_b_1*k_c_1*k_d_1*k_e_1*k_f_1*S_e_in_1;            %Marin Endurance Limit (gearset 2)

c = 1.34;
d = -0.085;

k_a_2 = c*(S_ut_1^d);                           %fatigue coefficients       
k_b_2 = 1;
k_c_2 = 1;
k_d_2 = 1;
k_e_2 = 1;
k_f_2 = 1.66;

S_e_2= k_a_2*k_b_2*k_c_2*k_d_2*k_e_2*k_f_2*S_e_in_2;       %Marin Endurance Limit (gearset 2)

sigmab_allowed_1 = S_e_1/n_s;                                 %Bending Stress (from marin endurance limit)
sigmab_allowed_2 = S_e_2/n_s;

T_1 = W_t_1*(d_p_1/24)*60;                                  %torques
T_2 = W_t_2*(d_p_2/24)*60;

PR_bending_1 = (T_1*omega_nc)/33000;                         %Gearset 1 & 2 power rating
PR_bending_2 = (T_2*omega_nc)/33000;

C_p_1 = 2100;
C_p_2 = 2300;

S_c_1 = (0.4*BH_1)-10;                                         %Surface Strength lim (kpsi)
S_c_2 = (0.4*BH_2)-10;

sigmac_allowed_1 = (-S_c_1/sqrt(n_s))*1000;                            %Surface Stress (kpsi)
sigmac_allowed_2 = (-S_c_2/sqrt(n_s))*1000;

r_1 = d_p_1/2;
r_2 = d_p_2/2;
r_3 = d_p_3/3;
r_4 = d_p_4/4;

W_t_c_1 = ((sigmac_allowed_1/-C_p_1)^2)*((r_1+r_2)*(F_1*cos(phi))/Kv_1)*57.3;   %surface transmitted loads (57.3 used to convert units)
W_t_c_2 = ((sigmac_allowed_2/-C_p_2)^2)*((r_3+r_4)*(F_3*cos(phi))/Kv_3)*57.3;

T_c_1 = W_t_c_1*r_1;
T_c_2 = W_t_c_2*r_3;

PR_surface_1 = (T_c_1*omega_nc)/33000;
PR_surface_2 = (T_c_2*omega_nc)/33000;

PRSET1 = min(PR_bending_1,PR_surface_1);                 %power rating set 1
PRSET2 = min(PR_bending_2,PR_surface_2);

PRGEARBOX = min(PRSET1,PRSET2);


fprintf('The power rating for the entire gear box is %s\n',PRGEARBOX)               %fancy
fprintf('The selected bending stress for gear set 1 is %s\n',sigmab_allowed_1)
fprintf('The selected bending stress for gear set 2 is %s\n',sigmab_allowed_2)
fprintf('The selected surface stress for gear set 1 is %s\n',sigmac_allowed_1)
fprintf('The selected surface stress for gear set 2 is %s\n',sigmac_allowed_2)






