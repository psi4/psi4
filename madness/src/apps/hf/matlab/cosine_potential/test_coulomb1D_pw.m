clear all;
close all;

alpha = 4.5;
npts = 30;
L = 30.0;

beta = L*L/(4*pi*pi);
gamma = 2*beta;

H = diag([-(npts-1)/2:(npts-1)/2].^2 - alpha*gamma*ones(1,npts));
H = H - alpha*beta*diag(ones(npts-1,1),-1) - alpha*beta*diag(ones(npts-1,1),1);
 
evH = eig(H)/gamma;
