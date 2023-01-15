function [E_out,C_out,VP_SCA,C_SCA,K_SCA,G_SCA] = run_SCADEM_joint(K,G,C,rho,phi_c,output)
% Rev. 2021-11-04 - Cosmetics
% Rev. 2021-11-01 - Important! Do not use version before.


% INPUT PARAMETERS
if isempty(output); output{1} = 'VP'; output{2} = 'R';  end

% START FRACTIONS
x1 = [(1-phi_c(1)) phi_c(1)];   x2 = [(1-phi_c(2)) phi_c(2)];   rho_b  = sum(x1.*rho);


% SELF-CONSISTENT APPROXIMATION (SCA)
[ E_opt,~,~,~] = opt_SCA_ul(K,G,C,x1);
[ ~,C_opt,~,~] = opt_SCA_ul(K,G,C,x2);
K_SCA  = E_opt(1);   
G_SCA  = E_opt(2); 
C_SCA  = C_opt;
VP_SCA = sqrt( (K_SCA + (4/3)*G_SCA)/rho_b);


% ELASTICS DEM
N_in = 1;               % First Inclusion (LHS of DEM)
[K_m1,G_m1,incl1] = dem(K_SCA,G_SCA,K(N_in),G(N_in),phi_c(1));
phi_DEM1   = flip(phi_c(1) - incl1'); 
K_DEM1     = flip(K_m1'); 
G_DEM1     = flip(G_m1');
f1         = [1-phi_c(1)+incl1'; flip(phi_DEM1)];
rho_b_DEM1 = flip(sum([f1(1,:).*rho(1); f1(2,:).*rho(2)],1)); 
VP_DEM1    = sqrt( (K_DEM1 + (4/3).*G_DEM1)./rho_b_DEM1);

N_in = 2;               % Second Inclusion (RHS of DEM)
[K_m2,G_m2,incl2] = dem(K_SCA,G_SCA,K(N_in),G(N_in),1-phi_c(1));
phi_DEM2   = phi_c(1) + incl2'; 
K_DEM2     = K_m2'; 
G_DEM2     = G_m2';
f2         = [1-phi_c(1)-incl2'; phi_DEM2+incl2];
rho_b_DEM2 = sum([f2(1,:).*rho(1); f2(2,:).*rho(2)],1); 
VP_DEM2    = sqrt( (K_DEM2 + (4/3).*G_DEM2)./rho_b_DEM2);


% ELASTICS OUTPUT
phi_DEM    = [phi_DEM1 phi_DEM2];
VP_DEM     = [VP_DEM1 VP_DEM2];
K_DEM      = [K_DEM1 K_DEM2];
G_DEM      = [G_DEM1 G_DEM2];

if strcmpi(output{1},'VP') == 1
    E_out  = [phi_DEM; VP_DEM]; E_out = E_out';    
    
elseif strcmpi(output{1},'Moduli') == 1
 	E_out  = [phi_DEM; K_DEM; G_DEM]; E_out = E_out';  
    
else
    error('Unkown Output Method for Elastics');
    
end


% ELECTRICS DEM
N_in = 1;               % First Inclusion (LHS of DEM) 
[C1,inclC1] = demC(C_SCA,C(N_in),phi_c(2));
C_phi_DEM1  = flip(phi_c(2) - inclC1'); C_DEM1 = flip(C1);

N_in = 2;               % Second Inclusion (RHS of DEM)
[C2,inclC2] = demC(C_SCA,C(N_in),1-phi_c(2));
C_phi_DEM2  = phi_c(2) + inclC2'; C_DEM2 = C2;


% ELECTRICS OUTPUT
C_phi_DEM   = [C_phi_DEM1 C_phi_DEM2];
C_DEM       = [C_DEM1' C_DEM2'];


if strcmpi(output{2},'R') == 1
    C_out   = [C_phi_DEM; 1./C_DEM]; C_out = C_out';
    
elseif strcmpi(output{2},'C') == 1
    C_out   = [C_phi_DEM; C_DEM]; C_out = C_out';
    
else
    error('Unkown Output Method for Electrics');
    
end

end

