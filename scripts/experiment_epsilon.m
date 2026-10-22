function experiment_epsilon()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
ratio_list = 0.5:0.5:5;
E_ratio = sweep_values(schemes, ratio_list, ...
    @(v) struct('eps_c', 1.5, 'eps_p', 1.5 / v), 1:n_seed, 'Epsilon ratio');
save(fullfile(repo, 'results', 'epsilon.mat'), 'schemes', 'ratio_list', 'E_ratio');
end
