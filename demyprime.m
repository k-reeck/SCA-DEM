function yprime=demyprime(t,y)
%function yprime=demyprime(t,y)
%used by DEM

%Written by T. Mukerji (2009), https://pangea.stanford.edu/research/srb/books/RPHII/RPHtools/
% Remarks: removed PQ calculation for arbitrary inclusions for speed
% reasons, can be found in original files 

global DEMINPT;

k1=DEMINPT(1); mu1=DEMINPT(2); k2=DEMINPT(3); mu2=DEMINPT(4); phic=DEMINPT(5); 

ka=k2; mua=mu2;
k=y(1); mu=y(2);
yprime=zeros(2,1);

% Modified for spherical only -> faster
CETA    = @(K,G) (G/6).*((9*K+8*G)./(K+2*G));
pa = (k + (4/3)*mu)./(ka + (4/3)*mu);
qa = (mu + CETA(k,mu))./(mua+CETA(k,mu));

krhs=(ka-k)*pa;
yprime(1)=krhs/(1-t);

murhs=(mua-mu)*qa;
yprime(2)=murhs/(1-t);
