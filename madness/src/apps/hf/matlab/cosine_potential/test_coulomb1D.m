clear all;
close all;

%% lattice constant
npts = 100;
L = 5.0;
x = linspace(-L/2, L/2, npts);
delx = x(2) - x(1);
%% 1-body potential
alpha = 12.5;
n = 1;
v = -alpha * (cos(2 * n * pi * x/ L) + 1);
V = diag(v);
%% kinetic energy
T = 2 * diag(ones(npts, 1)) - diag(ones(npts-1,1), 1) - diag(ones(npts-1,1), -1);
T(1,npts) = -1;
T(npts,1) = -1;
T = T ./ (delx)^2;
%% hamiltonian
H = T + V;