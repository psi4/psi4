clear all;
close all;

% lattice constant
npts = 15;
L = 30.0;
alpha = 13.5;

% build one component of the Hamiltonian
Hi = buildHamiltonian1D(L, npts,alpha);
ev1d = eig(Hi);

% build 3D Hamiltonian
x = linspace(-L/2, L/2, npts+1);
y = linspace(-L/2, L/2, npts+1);
z = linspace(-L/2, L/2, npts+1);
x = x(1:npts);
y = y(1:npts);
z = z(1:npts);
delx = x(2) - x(1);
% 1-body potential
v = zeros(npts,npts,npts);
for i = 1:npts
    for j = 1:npts
        for k = 1:npts
            v(i,j,k) = -alpha*(cos(2*pi*x(i)/ L) + cos(2*pi*y(j)/ L) + cos(2*pi*z(k)/L) + 3);
        end
    end
end
V = diag(v(:));
% kinetic energy
T = -0.5 * buildlaplacian3d(npts, delx);
% hamiltonian
H = T + V;
ev3d = eig(H);

error = abs(ev3d(1) - 3*ev1d(1));
fprintf(1,'ev1d(1) = %.8f  ev3d = %.8f  error = %.8f\n\n', ev1d(1), ev3d(1), error);