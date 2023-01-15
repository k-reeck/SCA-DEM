function [ E_opt,C_opt,E_res,C_res] = opt_SCA_ul(K,G,C,x)

E_opt = [];
C_opt = [];
E_res = [];
C_res = [];

% Hill's Average for Startmodel / by Archie for Conductivity
%K_m0 = 0.5*(sum(x.*K) + sum(x./K)^-1);
%G_m0 = 0.5*(sum(x.*G) + sum(x./G)^-1);

K_m0 = sum(x.*K);
G_m0 = sum(x.*G);
C_m0 = 0.5*(sum(x.*C) + sum(x./C)^-1);



% 1% Error estimate for 
%E_err = [K_m0/1e0 G_m0/1e0]; E0 = log([K_m0 G_m0]); 
E_err = [K_m0/1e2 G_m0/0.99e2]; E0 = log([K_m0 G_m0]); 

C_err = C_m0/1e2;
C0 = log(C_m0);

options1 = optimset('disp','off','FinDiffType','central','Algorithm','levenberg-marquardt');    
options2 = optimset('disp','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'MaxFunEvals',1000,'DiffMinChange',6e0,'DiffMaxChange',2e1,'Algorithm','levenberg-marquardt');      %Gute Einstellung
 
% Non-linear Solution
[E_opt,~,E_res] = lsqnonlin(@(E) opt_Elastics_lsqnonlin(E,E_err,K,G,x),E0,[],[],options1);     E_opt = exp(E_opt);
[C_opt,~,C_res] = lsqnonlin(@(C_m) opt_Cond_lsqnonlin(C_m,C_err,C,x),C0,[],[],options2);       C_opt = exp(C_opt);

% Linear Solution
%[E_opt, E_res]      = opt_Lin(K_m0,G_m0,K,G,x,asp,method);
%[C_opt, C_res,j]    = opt_Cond_Lin(C_m0,C,x,asp,method);

end

function e = opt_Elastics_lsqnonlin(par,err,K,G,x)
        
        K_m = exp(par(1));
        G_m = exp(par(2));
        

        CETA    = @(K,G) (G/6).*((9*K+8*G)./(K+2*G));
        P = (K_m + (4/3)*G_m)./(K + (4/3)*G_m);
        Q = (G_m + CETA(K_m,G_m))./(G+CETA(K_m,G_m));

        e1 = abs(sum(x.*(K - K_m).*P))./err(1);
        e2 = abs(sum(x.*(G - G_m).*Q))./err(2);
        
        e = sqrt(e1^2 + e2^2);
        %e = e1;
        %e = norm([sum(x.*(K - K_m).*P sum(x.*(G - G_m).*Q)]); 
        e1 = e;
                 
end

function e = opt_Cond_lsqnonlin(par,err,C,x)
        
        C_m = exp(par); 
        R = 1./(C+2*C_m);
        e = sum(x.*(C - C_m).*R);      
end

function [par_opt, res] = opt_Lin(K_m0,G_m0,K,G,x)
    
    e0 = 10;
    e  = K_m0;
    j  = 0;
    
    K_m = []; G_m = [];
    

    while (e > abs(e0)) && (j < 1000) 
    
        
        CETA    = @(K,G) (G/6).*((9*K+8*G)./(K+2*G));
        P = (K_m + (4/3)*G_m)./(K + (4/3)*G_m);
        Q = (G_m + CETA(K_m,G_m))./(G+CETA(K_m,G_m));
        
       
        K_m = sum(x.*K.*P)./sum(x.*P);   
        G_m = sum(x.*G.*Q)./sum(x.*Q);
        
        e    = norm([abs(K_m0 - K_m) abs(G_m0 - G_m)]);
        j    = j + 1;
        K_m0 = K_m;
        G_m0 = G_m;

    end
    
    if isempty(K_m)
        K_m = K_m0;
        G_m = G_m0;
    end
    
    par_opt = [K_m G_m];
    res = e;

end

function [par_opt,res,j] = opt_Cond_Lin(C_m0,C,x)
        
    e0 = C_m0/1e9;
    e  = C_m0;
    j = 0;
    C_m = []; 
    
    while (e > abs(e0)) && (j < 1000) 
        
       
        R = 1./(C+2*C_m);
        
        C_m = sum(x.*R)./sum(x.*C.*R);
        e    = abs(C_m0 - C_m);    	
        j    = j + 1;
        C_m0 = C_m;

    end
    
    
    if isempty(C_m)
        C_m = C_m0;       
        disp('C_m empty! Check data!');
    end
    
    par_opt = C_m;
    res = e;
end