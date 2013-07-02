clear all;
close all;

alpha = 12.5;
% npts has to be odd
npts = 51;
L = 5.0;

beta = L*L/(4*pi*pi);
gamma = 2*beta;
% range needs to go from -(npts-1)/2 ... (npts-1)/2
H = diag([-(npts-1)/2:(npts-1)/2].^2 - alpha*gamma*ones(1,npts));
H = H - alpha*beta*diag(ones(npts-1,1),-1) - alpha*beta*diag(ones(npts-1,1),1);
 
evH = eig(H)/gamma;
