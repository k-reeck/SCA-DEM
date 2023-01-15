function [C,por]=demC(C_m,Ci,phic)

% Modified by K. Reeck (kreeck@geomar.de) for electric SCA-DEM formulation
% for T. Mukerji (2009), RPHtools Matlab Scripts. https://pangea.stanford.edu/research/srb/books/RPHII/RPHtools/

global DEMINPT;
DEMINPT=ones(1,4);
DEMINPT(1)=C_m; 
DEMINPT(2)=Ci; 
DEMINPT(3)=phic; 

[tout, yout]=ode45m('demCyprime',0.00,0.99999,C_m,1e-10);

 C=real(yout(:,1)); 
por=phic*tout;

clear DEMINPT;
