function call_SCADEM_test()

% Replicate Clay/Water Mixture for phi_c = [0.5 0.5] according to 
% Figure 2b in Manuscript, inclusions AR are fixed to spherical for faster
% computing

% Bulk Modulus
K.C     = 20.9*1e9; 
K.W     = 2.29*1e9; 

% Shear Modulus
G.C     = 6.85.*1e9;         
G.W     = 1e-5;

% Conductivity
C.C     = 0.02; 
C.W     = 3.176;            

% Density
rho.C   = 2.58*1e3; 
rho.W   = 1.025*1e3;


    output = {'Moduli' 'C'};
    
    [E_out,C_out,~,~,~,~] = run_SCADEM_joint([K.C K.W],[G.C G.W],[C.C C.W],[rho.C rho.W],[0.5 0.5],output);

    

    
    figure;
    yyaxis left
    plot(E_out(:,1).*100,E_out(:,2)./1e9,'b-',E_out(:,1).*100,E_out(:,3)./1e9,'b--'); grid on;
    xlabel('Concentration x_2 (%)');
    ylabel('Elastic Moduli (GPa)');
    yyaxis right
    semilogy(C_out(:,1).*100,C_out(:,2),'r-'); 
    ylabel('Electric Conductivity (S/m)');
    legend('Bulk Modulus','Shear Modulus','Conductivity');
    



end

