function experiment_ap_number()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
M_list = 3:12;
E_M = zeros(numel(schemes), numel(M_list));
rate_M = zeros(numel(schemes), numel(M_list));
for s = 1:n_seed
    for j = 1:numel(M_list)
        fprintf('APs, M = %d, seed %d\n', M_list(j), s);
        sc = build_scenario(M_list(j), 5, s);
        for i = 1:numel(schemes)
            metrics = optimize_ucmec(sc, schemes{i}).metrics;
            E_M(i, j) = E_M(i, j) + metrics.E_tot;
            rate_M(i, j) = rate_M(i, j) + metrics.rate_Mbps;
        end
    end
end
E_M = E_M / n_seed;
rate_M = rate_M / n_seed;
save(fullfile(repo, 'results', 'ap_number.mat'), 'schemes', 'M_list', 'E_M', 'rate_M');
end
