function experiment_user_number()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
N_list = 5:14;
n_user = numel(N_list);
n_scheme = numel(schemes);
E_off = zeros(n_scheme, n_user);
E_con = zeros(n_scheme, n_user);
E_tot = zeros(n_scheme, n_user);
D_off = zeros(n_scheme, n_user);
D_con = zeros(n_scheme, n_user);
D_tot = zeros(n_scheme, n_user);
rate = zeros(n_scheme, n_user);

for s = 1:n_seed
    for j = 1:n_user
        fprintf('Users, N = %d, seed %d\n', N_list(j), s);
        sc = build_scenario(5, N_list(j), s);
        for i = 1:n_scheme
            m = optimize_ucmec(sc, schemes{i}).metrics;
            E_off(i, j) = E_off(i, j) + m.E_off;
            E_con(i, j) = E_con(i, j) + m.E_con;
            E_tot(i, j) = E_tot(i, j) + m.E_tot;
            D_off(i, j) = D_off(i, j) + m.D_off_ms;
            D_con(i, j) = D_con(i, j) + m.D_con_ms;
            D_tot(i, j) = D_tot(i, j) + m.D_tot_ms;
            rate(i, j) = rate(i, j) + m.rate_Mbps;
        end
    end
end
E_off = E_off / n_seed;
E_con = E_con / n_seed;
E_tot = E_tot / n_seed;
D_off = D_off / n_seed;
D_con = D_con / n_seed;
D_tot = D_tot / n_seed;
rate = rate / n_seed;

save(fullfile(repo, 'results', 'user_number.mat'), ...
    'N_list', 'schemes', 'E_off', 'E_con', 'E_tot', ...
    'D_off', 'D_con', 'D_tot', 'rate');
end
