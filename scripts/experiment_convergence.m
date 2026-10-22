function experiment_convergence()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, ~, n_seed] = setup_path();
out_dir = fullfile(repo, 'results');
M_list = 3:6;
q_list = [10, 50, 100, 200];
tmax = 40;
seeds = 1:n_seed;

convergence_ap = nan(numel(M_list), tmax);
for i = 1:numel(M_list)
    fprintf('Convergence, M = %d\n', M_list(i));
    hist = nan(n_seed, tmax);
    for s = 1:n_seed
        sc = build_scenario(M_list(i), 5, seeds(s));
        sc.tmax = tmax;
        y = optimize_ucmec(sc, 'proposed').history.E_tot;
        y(end + 1:tmax) = y(end);
        hist(s, :) = y(1:tmax);
    end
    convergence_ap(i, :) = mean(hist, 1);
end
save(fullfile(out_dir, 'convergence_ap.mat'), 'convergence_ap', 'M_list', 'tmax');

convergence_q = nan(numel(q_list), tmax);
for i = 1:numel(q_list)
    fprintf('Convergence, q = %g\n', q_list(i));
    hist = nan(n_seed, tmax);
    for s = 1:n_seed
        sc = build_scenario(5, 5, seeds(s), struct('q', q_list(i)));
        sc.tmax = tmax;
        y = optimize_ucmec(sc, 'proposed').history.E_tot;
        y(end + 1:tmax) = y(end);
        hist(s, :) = y(1:tmax);
    end
    convergence_q(i, :) = mean(hist, 1);
end
save(fullfile(out_dir, 'convergence_q.mat'), 'convergence_q', 'q_list', 'tmax');
end
