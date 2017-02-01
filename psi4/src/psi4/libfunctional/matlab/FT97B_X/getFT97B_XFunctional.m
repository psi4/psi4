%
% @BEGIN LICENSE
%
% Psi4: an open-source quantum chemistry software package
%
% Copyright (c) 2007-2017 The Psi4 Developers.
%
% The copyrights for code used from other parties are included in
% the corresponding files.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% @END LICENSE
%

function data = getFT97B_Functional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'd0', 'd1', 'd2'};
data.param_vals = [-3/8*3^(1/3)*4^(2/3)*pi^(-1/3), 0.002913644, 0.0009474169, 2501.149];

syms c d0 d1 d2 real;
syms chi2 n s;
d = d0 + d1*s/(d2^2 + s);
ldax = c*n^(4/3);
FX = 1-d*chi2/(c*sqrt(1+9*d^2*chi2*asinh(chi2)^2));

data.functional = subs(ldax,n,rho_a) * subs(FX, {chi2, n,s}, {gamma_aa/rho_a^(8/3),rho_a,gamma_aa}) + ...
                  subs(ldax,n,rho_b) * subs(FX, {chi2, n,s}, {gamma_bb/rho_b^(8/3),rho_b,gamma_bb}); 
data.functional_a0 = subs(ldax,n,rho_b) * subs(FX, {chi2, n,s}, {gamma_bb/rho_b^(8/3),rho_b,gamma_bb});  
data.functional_b0 = subs(ldax,n,rho_a) * subs(FX, {chi2, n,s}, {gamma_aa/rho_a^(8/3),rho_a,gamma_aa}); 
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 1;

data.name = 'FT97B_X';
data.citation = 'M. Filatov and W. Theil, Mol. Phys., 91(5), 847-859, 1997.';
data.description = 'Filitov and Theil 1997 Exchange';
