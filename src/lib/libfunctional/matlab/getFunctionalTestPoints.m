function points = getFunctionalTestPoints()

points.n = 23;

rho_a_ = zeros(points.n,1);
rho_b_ = zeros(points.n,1);
gamma_aa_ = zeros(points.n,1);
gamma_ab_ = zeros(points.n,1);
gamma_bb_ = zeros(points.n,1);
tau_a_ = zeros(points.n,1);
tau_b_ = zeros(points.n,1);

index = 1;
rho_a_(index) = 0.17E+01; rho_b_(index) = 0.17E+01; gamma_aa_(index) = 0.81E-11; gamma_ab_(index) = 0.81E-11; gamma_bb_(index) = 0.81E-11; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.17E+01; rho_b_(index) = 0.17E+01; gamma_aa_(index) = 0.17E+01; gamma_ab_(index) = 0.17E+01; gamma_bb_(index) = 0.17E+01; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.15E+01; rho_b_(index) = 0.15E+01; gamma_aa_(index) = 0.36E+02; gamma_ab_(index) = 0.36E+02; gamma_bb_(index) = 0.36E+02; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.88E-01; rho_b_(index) = 0.88E-01; gamma_aa_(index) = 0.87E-01; gamma_ab_(index) = 0.87E-01; gamma_bb_(index) = 0.87E-01; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.18E+04; rho_b_(index) = 0.18E+04; gamma_aa_(index) = 0.55E+00; gamma_ab_(index) = 0.55E+00; gamma_bb_(index) = 0.55E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.18E+04; rho_b_(index) = 0.18E+04; gamma_aa_(index) = 0.86E+04; gamma_ab_(index) = 0.86E+04; gamma_bb_(index) = 0.86E+04; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.16E+04; rho_b_(index) = 0.16E+04; gamma_aa_(index) = 0.37E+10; gamma_ab_(index) = 0.37E+10; gamma_bb_(index) = 0.37E+10; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.26E+00; rho_b_(index) = 0.26E+00; gamma_aa_(index) = 0.28E+00; gamma_ab_(index) = 0.28E+00; gamma_bb_(index) = 0.28E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.53E+05; rho_b_(index) = 0.53E+05; gamma_aa_(index) = 0.96E+05; gamma_ab_(index) = 0.96E+05; gamma_bb_(index) = 0.96E+05; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.47E+05; rho_b_(index) = 0.47E+05; gamma_aa_(index) = 0.29E+14; gamma_ab_(index) = 0.29E+14; gamma_bb_(index) = 0.29E+14; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.15E+00; rho_b_(index) = 0.15E+00; gamma_aa_(index) = 0.16E+00; gamma_ab_(index) = 0.16E+00; gamma_bb_(index) = 0.16E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.35E+01; rho_b_(index) = 0.00E+00; gamma_aa_(index) = 0.46E-10; gamma_ab_(index) = 0.00E+00; gamma_bb_(index) = 0.00E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.35E+01; rho_b_(index) = 0.00E+00; gamma_aa_(index) = 0.34E+01; gamma_ab_(index) = 0.00E+00; gamma_bb_(index) = 0.00E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.30E+01; rho_b_(index) = 0.00E+00; gamma_aa_(index) = 0.20E+03; gamma_ab_(index) = 0.00E+00; gamma_bb_(index) = 0.00E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.58E-01; rho_b_(index) = 0.00E+00; gamma_aa_(index) = 0.47E-01; gamma_ab_(index) = 0.00E+00; gamma_bb_(index) = 0.00E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.82E+02; rho_b_(index) = 0.81E+02; gamma_aa_(index) = 0.49E+07; gamma_ab_(index) = 0.49E+07; gamma_bb_(index) = 0.49E+07; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.39E+02; rho_b_(index) = 0.38E+02; gamma_aa_(index) = 0.81E+06; gamma_ab_(index) = 0.82E+06; gamma_bb_(index) = 0.82E+06; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.13E+00; rho_b_(index) = 0.95E-01; gamma_aa_(index) = 0.15E+00; gamma_ab_(index) = 0.18E+00; gamma_bb_(index) = 0.22E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.78E-01; rho_b_(index) = 0.31E-01; gamma_aa_(index) = 0.41E-02; gamma_ab_(index) = 0.38E-02; gamma_bb_(index) = 0.36E-02; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.50E+02; rho_b_(index) = 0.49E+02; gamma_aa_(index) = 0.11E+06; gamma_ab_(index) = 0.11E+06; gamma_bb_(index) = 0.11E+06; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.40E+02; rho_b_(index) = 0.40E+02; gamma_aa_(index) = 0.99E+05; gamma_ab_(index) = 0.98E+05; gamma_bb_(index) = 0.98E+05; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.12E+00; rho_b_(index) = 0.10E+00; gamma_aa_(index) = 0.12E+00; gamma_ab_(index) = 0.13E+00; gamma_bb_(index) = 0.14E+00; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;
    rho_a_(index) = 0.48E-01; rho_b_(index) = 0.25E-01; gamma_aa_(index) = 0.46E-02; gamma_ab_(index) = 0.44E-02; gamma_bb_(index) = 0.41E-02; tau_a_(index) = 0.00E+00; tau_b_(index) = 0.00E+00; index = index + 1;

points.rho_a = rho_a_;
points.rho_b = rho_b_;
points.gamma_aa = gamma_aa_;
points.gamma_ab = gamma_ab_;
points.gamma_bb = gamma_bb_;
points.tau_a = tau_a_;
points.tau_b = tau_b_;

