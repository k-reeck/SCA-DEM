function yprime=demCyprime(t,y)
% Modified by K. Reeck (kreeck@geomar.de) for electric SCA-DEM formulation
% for T. Mukerji (2009), RPHtools Matlab Scripts. https://pangea.stanford.edu/research/srb/books/RPHII/RPHtools/
%
% Remarks: Removed R-calculation for arbitrary inclusiosn for speed reasons
% can be changed by including other demag factors (L) after Osborne (1945). 

global DEMINPT;

C_m  = DEMINPT(1); 
Ci   = DEMINPT(2);  
phic =DEMINPT(3);


C = y(1);
yprime=zeros(1,1);


L = [1/3; 1/3; 1/3];
Ri = (1/9)*sum( 1./(L*Ci + (1-L)*C) );



yprime(1) = ((Ci-C)*(3*C)*Ri)/(1-t);


