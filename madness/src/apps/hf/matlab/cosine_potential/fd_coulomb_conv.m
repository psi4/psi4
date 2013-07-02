% This script tests a finite difference scheme on the 3D coulomb
% convolution problem.

clear all;
close all;

% Size of cube
nsize = 20;
% Dimensions of cube
L = 5.0;
% points in real space
x = linspace(-L/2,L/2,nsize);
y = linspace(-L/2,L/2,nsize);
z = linspace(-L/2,L/2,nsize);
% make rho
fprintf(1,'making rho ...\n\n');
rho = zeros(nsize,nsize,nsize);
for ri = 1:nsize
    for rj = 1:nsize
        for rk = 1:nsize
            rho(ri,rj,rk) = exp(-5*(x(ri)*x(ri) + y(rj)*y(rj) + z(rk)*z(rk)));
        end
    end
end
% compute indicies
[I,J,K] = ind2sub([nsize nsize nsize], [1:nsize*nsize*nsize]);

% return value
rv = zeros(nsize,nsize,nsize);

% main loop
% if you think of this convolution as a matrix-vector multiply,
% mi and mj index into the matrix which has the size of [npoints, npoints]
% where npoints is the total number of points in the cube
fprintf(1,'running main loop ...\n\n');
for mi = 1:nsize*nsize*nsize
    fprintf(1,'mi = %d\n', mi);
    for mj = 1:nsize*nsize*nsize
        % compute indices to express mi and mj into a 3-D index
        % so we can associate ki and kj with variables x and x'
        ki = [I(mi) J(mi) K(mi)]; 
        kj = [I(mj) J(mj) K(mj)]; 
        % compute distance
        distx = x(ki(1)) - x(kj(1));
        disty = x(ki(2)) - x(kj(2));
        distz = x(ki(3)) - x(kj(3));
        dist = sqrt(distx^2 + disty^2 + distz^2);
        temp = rho(kj(1), kj(2), kj(3))/dist;
        rv(ki) = rv(ki) + temp;
    end
end
