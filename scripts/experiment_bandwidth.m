function experiment_bandwidth()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
B_list = 10:10:50;
E_B = sweep_values(schemes, B_list, @(v) struct('B_MHz', v), 1:n_seed, 'Bandwidth');
save(fullfile(repo, 'results', 'bandwidth.mat'), 'schemes', 'B_list', 'E_B');
end
