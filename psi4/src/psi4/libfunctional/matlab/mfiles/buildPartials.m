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

function buildPartials(fun, name, deriv)

syms rho_a rho_b gamma_aa gamma_ab gamma_bb tau_a tau_b real;
zero = sym('0');

if (deriv >=0)
    ccode(fun + zero, 'file', 'v');
end

if (deriv >= 1)
    ccode(diff(fun + zero,rho_a,1),'file','v_rho_a');
    ccode(diff(fun + zero,rho_b,1),'file','v_rho_b');
    ccode(diff(fun + zero,gamma_aa,1),'file','v_gamma_aa');
    ccode(diff(fun + zero,gamma_ab,1),'file','v_gamma_ab');
    ccode(diff(fun + zero,gamma_bb,1),'file','v_gamma_bb');
    ccode(diff(fun + zero,tau_a,1),'file','v_tau_a');
    ccode(diff(fun + zero,tau_b,1),'file','v_tau_b');
end

if (deriv >= 2)
    ccode(diff(diff(fun + zero,rho_a,1),rho_a,1),'file','v_rho_a_rho_a');
    ccode(diff(diff(fun + zero,rho_a,1),rho_b,1),'file','v_rho_a_rho_b');
    ccode(diff(diff(fun + zero,rho_b,1),rho_b,1),'file','v_rho_b_rho_b');
    ccode(diff(diff(fun + zero,rho_a,1),gamma_aa,1),'file','v_rho_a_gamma_aa');
    ccode(diff(diff(fun + zero,rho_a,1),gamma_ab,1),'file','v_rho_a_gamma_ab');
    ccode(diff(diff(fun + zero,rho_a,1),gamma_bb,1),'file','v_rho_a_gamma_bb');
    ccode(diff(diff(fun + zero,rho_b,1),gamma_aa,1),'file','v_rho_b_gamma_aa');
    ccode(diff(diff(fun + zero,rho_b,1),gamma_ab,1),'file','v_rho_b_gamma_ab');
    ccode(diff(diff(fun + zero,rho_b,1),gamma_bb,1),'file','v_rho_b_gamma_bb');
    ccode(diff(diff(fun + zero,gamma_aa,1),gamma_aa,1),'file','v_gamma_aa_gamma_aa');
    ccode(diff(diff(fun + zero,gamma_aa,1),gamma_ab,1),'file','v_gamma_aa_gamma_ab');
    ccode(diff(diff(fun + zero,gamma_aa,1),gamma_bb,1),'file','v_gamma_aa_gamma_bb');
    ccode(diff(diff(fun + zero,gamma_ab,1),gamma_ab,1),'file','v_gamma_ab_gamma_ab');
    ccode(diff(diff(fun + zero,gamma_ab,1),gamma_bb,1),'file','v_gamma_ab_gamma_bb');
    ccode(diff(diff(fun + zero,gamma_bb,1),gamma_bb,1),'file','v_gamma_bb_gamma_bb');
    ccode(diff(diff(fun + zero,rho_a,1),tau_a,1),'file','v_rho_a_tau_a');
    ccode(diff(diff(fun + zero,rho_a,1),tau_b,1),'file','v_rho_a_tau_b');
    ccode(diff(diff(fun + zero,rho_b,1),tau_a,1),'file','v_rho_b_tau_a');
    ccode(diff(diff(fun + zero,rho_b,1),tau_b,1),'file','v_rho_b_tau_b');
    ccode(diff(diff(fun + zero,tau_a,1),tau_a,1),'file','v_tau_a_tau_a');
    ccode(diff(diff(fun + zero,tau_a,1),tau_b,1),'file','v_tau_a_tau_b');
    ccode(diff(diff(fun + zero,tau_b,1),tau_b,1),'file','v_tau_b_tau_b');
    ccode(diff(diff(fun + zero,gamma_aa,1),tau_a,1),'file','v_gamma_aa_tau_a');
    ccode(diff(diff(fun + zero,gamma_ab,1),tau_a,1),'file','v_gamma_ab_tau_a');
    ccode(diff(diff(fun + zero,gamma_bb,1),tau_a,1),'file','v_gamma_bb_tau_a');
    ccode(diff(diff(fun + zero,gamma_aa,1),tau_b,1),'file','v_gamma_aa_tau_b');
    ccode(diff(diff(fun + zero,gamma_ab,1),tau_b,1),'file','v_gamma_ab_tau_b');
    ccode(diff(diff(fun + zero,gamma_bb,1),tau_b,1),'file','v_gamma_bb_tau_b');
end

system(['./purify.py ' name]);
system(['./expei.py ' name]);
