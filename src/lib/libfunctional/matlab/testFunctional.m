function testFunctional(data, weights,outfile)

out = fopen(outfile,'w');
tol = 1E-20;

fprintf(out,'Testing %s Functional\n\n', data.name);

fprintf(out,' Description: \n');
fprintf(out,'---------------\n');
fprintf(out,'%s\n\n',data.description);

fprintf(out,' Citation: \n');
fprintf(out,'---------------\n');
fprintf(out,'%s\n\n',data.citation);

functional = weights(1)*data(1).functional;
functional_a0 = weights(1)*data(1).functional_a0;
functional_b0 = weights(1)*data(1).functional_b0;
functional_a0b0 = weights(1)*data(1).functional_a0b0;
for k = 2:length(data)
    functional = functional + weights(k)*data(k).functional;
    functional_a0 = functional_a0 + weights(k)*data(k).functional_a0;
    functional_b0 = functional_b0 + weights(k)*data(k).functional_b0;
    functional_a0b0 = functional_a0b0 + weights(k)*data(k).functional_a0b0;
end

for l = 1:length(data)
    for k = 1:length(data(l).param_names)
        functional = subs(functional,sym(data(l).param_names{k}),data(l).param_vals(k));
        functional_a0 = subs(functional_a0,sym(data(l).param_names{k}),data(l).param_vals(k));
        functional_b0 = subs(functional_b0,sym(data(l).param_names{k}),data(l).param_vals(k));
        functional_a0b0 = subs(functional_a0b0,sym(data(l).param_names{k}),data(l).param_vals(k));
    end
end

points = getFunctionalTestPoints();
syms rho_a rho_b gamma_aa gamma_ab gamma_bb tau_a tau_b real;


for k = 1:points.n
    ra = points.rho_a(k);
    rb = points.rho_b(k);
    gaa = points.gamma_aa(k);
    gab = points.gamma_ab(k);
    gbb = points.gamma_bb(k);
    ta = points.tau_a(k);
    tb = points.tau_b(k);

    fprintf(out, ' **rho_a= %8.2E rho_b= %8.2E gamma_aa= %8.2E gamma_ab= %8.2E gamma_bb= %8.2E tau_a= %8.2E tau_b= %8.2E\n\n',...
        ra, rb, gaa, gab, gbb, ta, tb);        

    
    if (ra > tol && rb > tol)
        Q = functional;
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','functional',double(val));

        fprintf(out,'\n');       
 
        Q = diff(functional,rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a',double(val));

        Q = diff(functional,rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b',double(val));

        Q = diff(functional,gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa',double(val));

        Q = diff(functional,gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab',double(val));

        Q = diff(functional,gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb',double(val));

        Q = diff(functional,tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a',double(val));

        Q = diff(functional,tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b',double(val));
        
        fprintf(out,'\n');       

        Q = diff(diff(functional,rho_a),rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_a',double(val));

        Q = diff(diff(functional,rho_a),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_b',double(val));

        Q = diff(diff(functional,rho_b),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_rho_b',double(val));

        Q = diff(diff(functional,rho_a),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_aa',double(val));

        Q = diff(diff(functional,rho_a),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_ab',double(val));

        Q = diff(diff(functional,rho_a),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_bb',double(val));

        Q = diff(diff(functional,rho_b),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_aa',double(val));

        Q = diff(diff(functional,rho_b),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_ab',double(val));

        Q = diff(diff(functional,rho_b),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_bb',double(val));

        Q = diff(diff(functional,gamma_aa),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_aa',double(val));

        Q = diff(diff(functional,gamma_aa),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_ab',double(val));

        Q = diff(diff(functional,gamma_aa),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_bb',double(val));

        Q = diff(diff(functional,gamma_ab),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_ab',double(val));

        Q = diff(diff(functional,gamma_ab),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_bb',double(val));

        Q = diff(diff(functional,gamma_bb),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_gamma_bb',double(val));

        Q = diff(diff(functional,rho_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_a',double(val));

        Q = diff(diff(functional,rho_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_b',double(val));

        Q = diff(diff(functional,rho_b),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_a',double(val));

        Q = diff(diff(functional,rho_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_b',double(val));

        Q = diff(diff(functional,tau_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_a',double(val));

        Q = diff(diff(functional,tau_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_b',double(val));

        Q = diff(diff(functional,tau_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b_tau_b',double(val));

        Q = diff(diff(functional,gamma_aa),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_a',double(val));

        Q = diff(diff(functional,gamma_aa),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_b',double(val));

        Q = diff(diff(functional,gamma_ab),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_a',double(val));

        Q = diff(diff(functional,gamma_ab),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_b',double(val));

        Q = diff(diff(functional,gamma_bb),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_a',double(val));

        Q = diff(diff(functional,gamma_bb),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_b',double(val));

    elseif (ra > tol) 
        Q = functional_b0;
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','functional',double(val));

        fprintf(out,'\n');       
 
        Q = diff(functional_b0,rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a',double(val));

        Q = diff(functional_b0,rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b',double(val));

        Q = diff(functional_b0,gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa',double(val));

        Q = diff(functional_b0,gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab',double(val));

        Q = diff(functional_b0,gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb',double(val));

        Q = diff(functional_b0,tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a',double(val));

        Q = diff(functional_b0,tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b',double(val));
        
        fprintf(out,'\n');       

        Q = diff(diff(functional_b0,rho_a),rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_a',double(val));

        Q = diff(diff(functional_b0,rho_a),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_b',double(val));

        Q = diff(diff(functional_b0,rho_b),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_rho_b',double(val));

        Q = diff(diff(functional_b0,rho_a),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_aa',double(val));

        Q = diff(diff(functional_b0,rho_a),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_ab',double(val));

        Q = diff(diff(functional_b0,rho_a),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_bb',double(val));

        Q = diff(diff(functional_b0,rho_b),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_aa',double(val));

        Q = diff(diff(functional_b0,rho_b),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_ab',double(val));

        Q = diff(diff(functional_b0,rho_b),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_bb',double(val));

        Q = diff(diff(functional_b0,gamma_aa),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_aa',double(val));

        Q = diff(diff(functional_b0,gamma_aa),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_ab',double(val));

        Q = diff(diff(functional_b0,gamma_aa),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_bb',double(val));

        Q = diff(diff(functional_b0,gamma_ab),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_ab',double(val));

        Q = diff(diff(functional_b0,gamma_ab),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_bb',double(val));

        Q = diff(diff(functional_b0,gamma_bb),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_gamma_bb',double(val));

        Q = diff(diff(functional_b0,rho_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_a',double(val));

        Q = diff(diff(functional_b0,rho_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_b',double(val));

        Q = diff(diff(functional_b0,rho_b),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_a',double(val));

        Q = diff(diff(functional_b0,rho_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_b',double(val));

        Q = diff(diff(functional_b0,tau_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_a',double(val));

        Q = diff(diff(functional_b0,tau_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_b',double(val));

        Q = diff(diff(functional_b0,tau_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b_tau_b',double(val));

        Q = diff(diff(functional_b0,gamma_aa),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_a',double(val));

        Q = diff(diff(functional_b0,gamma_aa),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_b',double(val));

        Q = diff(diff(functional_b0,gamma_ab),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_a',double(val));

        Q = diff(diff(functional_b0,gamma_ab),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_b',double(val));

        Q = diff(diff(functional_b0,gamma_bb),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_a',double(val));

        Q = diff(diff(functional_b0,gamma_bb),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_b',double(val));

    elseif (rb > tol) 
        Q = functional_a0;
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','functional',double(val));

        fprintf(out,'\n');       
 
        Q = diff(functional_a0,rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a',double(val));

        Q = diff(functional_a0,rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b',double(val));

        Q = diff(functional_a0,gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa',double(val));

        Q = diff(functional_a0,gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab',double(val));

        Q = diff(functional_a0,gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb',double(val));

        Q = diff(functional_a0,tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a',double(val));

        Q = diff(functional_a0,tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b',double(val));
        
        fprintf(out,'\n');       

        Q = diff(diff(functional_a0,rho_a),rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_a',double(val));

        Q = diff(diff(functional_a0,rho_a),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_b',double(val));

        Q = diff(diff(functional_a0,rho_b),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_rho_b',double(val));

        Q = diff(diff(functional_a0,rho_a),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_aa',double(val));

        Q = diff(diff(functional_a0,rho_a),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_ab',double(val));

        Q = diff(diff(functional_a0,rho_a),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_bb',double(val));

        Q = diff(diff(functional_a0,rho_b),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_aa',double(val));

        Q = diff(diff(functional_a0,rho_b),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_ab',double(val));

        Q = diff(diff(functional_a0,rho_b),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_bb',double(val));

        Q = diff(diff(functional_a0,gamma_aa),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_aa',double(val));

        Q = diff(diff(functional_a0,gamma_aa),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_ab',double(val));

        Q = diff(diff(functional_a0,gamma_ab),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_ab',double(val));

        Q = diff(diff(functional_a0,gamma_ab),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_bb',double(val));

        Q = diff(diff(functional_a0,gamma_bb),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_gamma_bb',double(val));

        Q = diff(diff(functional_a0,gamma_aa),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_bb',double(val));

        Q = diff(diff(functional_a0,rho_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_a',double(val));

        Q = diff(diff(functional_a0,rho_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_b',double(val));

        Q = diff(diff(functional_a0,rho_b),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_a',double(val));

        Q = diff(diff(functional_a0,rho_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_b',double(val));

        Q = diff(diff(functional_a0,tau_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_a',double(val));

        Q = diff(diff(functional_a0,tau_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_b',double(val));

        Q = diff(diff(functional_a0,tau_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b_tau_b',double(val));

        Q = diff(diff(functional_a0,gamma_aa),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_a',double(val));

        Q = diff(diff(functional_a0,gamma_aa),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_b',double(val));

        Q = diff(diff(functional_a0,gamma_ab),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_a',double(val));

        Q = diff(diff(functional_a0,gamma_ab),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_b',double(val));

        Q = diff(diff(functional_a0,gamma_bb),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_a',double(val));

        Q = diff(diff(functional_a0,gamma_bb),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_b',double(val));

    else  
        Q = functional_a0b0;
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','functional',double(val));

        fprintf(out,'\n');       
 
        Q = diff(functional_a0b0,rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a',double(val));

        Q = diff(functional_a0b0,rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b',double(val));

        Q = diff(functional_a0b0,gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa',double(val));

        Q = diff(functional_a0b0,gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab',double(val));

        Q = diff(functional_a0b0,gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb',double(val));

        Q = diff(functional_a0b0,tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a',double(val));

        Q = diff(functional_a0b0,tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b',double(val));
        
        fprintf(out,'\n');       

        Q = diff(diff(functional_a0b0,rho_a),rho_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_a',double(val));

        Q = diff(diff(functional_a0b0,rho_a),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_rho_b',double(val));

        Q = diff(diff(functional_a0b0,rho_b),rho_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_rho_b',double(val));

        Q = diff(diff(functional_a0b0,rho_a),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_aa',double(val));

        Q = diff(diff(functional_a0b0,rho_a),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_ab',double(val));

        Q = diff(diff(functional_a0b0,rho_a),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_gamma_bb',double(val));

        Q = diff(diff(functional_a0b0,rho_b),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_aa',double(val));

        Q = diff(diff(functional_a0b0,rho_b),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_ab',double(val));

        Q = diff(diff(functional_a0b0,rho_b),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_gamma_bb',double(val));

        Q = diff(diff(functional_a0b0,gamma_aa),gamma_aa);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_aa',double(val));

        Q = diff(diff(functional_a0b0,gamma_aa),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_ab',double(val));

        Q = diff(diff(functional_a0b0,gamma_aa),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_gamma_bb',double(val));

        Q = diff(diff(functional_a0b0,gamma_ab),gamma_ab);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_ab',double(val));

        Q = diff(diff(functional_a0b0,gamma_ab),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_gamma_bb',double(val));

        Q = diff(diff(functional_a0b0,gamma_bb),gamma_bb);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_gamma_bb',double(val));

        Q = diff(diff(functional_a0b0,rho_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_a',double(val));

        Q = diff(diff(functional_a0b0,rho_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_a_tau_b',double(val));

        Q = diff(diff(functional_a0b0,rho_b),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_a',double(val));

        Q = diff(diff(functional_a0b0,rho_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_rho_b_tau_b',double(val));

        Q = diff(diff(functional_a0b0,tau_a),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_a',double(val));

        Q = diff(diff(functional_a0b0,tau_a),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_a_tau_b',double(val));

        Q = diff(diff(functional_a0b0,tau_b),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_tau_b_tau_b',double(val));

        Q = diff(diff(functional_a0b0,gamma_aa),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_a',double(val));

        Q = diff(diff(functional_a0b0,gamma_aa),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_aa_tau_b',double(val));

        Q = diff(diff(functional_a0b0,gamma_ab),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_a',double(val));

        Q = diff(diff(functional_a0b0,gamma_ab),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_ab_tau_b',double(val));

        Q = diff(diff(functional_a0b0,gamma_bb),tau_a);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_a',double(val));

        Q = diff(diff(functional_a0b0,gamma_bb),tau_b);
        val = subs(Q,[rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b], ...
            [ra rb gaa gab gbb ta tb]);
        fprintf(out,'  %-20s= %30.11E\n','v_gamma_bb_tau_b',double(val));

    end 
    fprintf(out,'\n');
end

fclose(out);
