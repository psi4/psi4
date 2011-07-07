function H = buildHamiltonian1D(L, npts, alpha)

% lattice constant
x = linspace(-L/2, L/2, npts+1);
x = x(1:npts);
delx = x(2) - x(1);
% 1-body potential
n = 1;
v = -alpha * (cos(2 * n * pi * x/ L) + 1);
V = diag(v);
% kinetic energy
T = diag(ones(npts, 1)) - 0.5*diag(ones(npts-1,1), 1) - 0.5*diag(ones(npts-1,1), -1);
T(1,npts) = -0.5;
T(npts,1) = -0.5;
T = T ./ (delx^2);
% hamiltonian
H = T + V;