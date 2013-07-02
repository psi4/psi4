clear all;
close all;

%% lattice constant
npts = 10;
L = 5.0;
x = linspace(-L/2, L/2, npts+1);
y = linspace(-L/2, L/2, npts+1);
z = linspace(-L/2, L/2, npts+1);
x = x(1:npts);
y = y(1:npts);
z = z(1:npts);
delx = x(2) - x(1);
%% 1-body potential
alpha = 12.5;
v = zeros(npts,npts,npts);
for i = 1:npts
    for j = 1:npts
        for k = 1:npts
            v(i,j,k) = -alpha*(cos(2*pi*x(i)/ L)*cos(2*pi*y(j)/ L)*cos(2*pi*z(k)/L) + 1);
        end
    end
end
V = diag(v(:));
%% kinetic energy
T = -0.5 * buildlaplacian3d(npts, delx);
%% hamiltonian
H = T + V;