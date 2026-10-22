function experiment_cpu()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
C_list = 10:10:50;
E_C = sweep_values(schemes, C_list, @(v) struct('C_GHz', v), 1:n_seed, 'CPU');
save(fullfile(repo, 'results', 'cpu.mat'), 'schemes', 'C_list', 'E_C');
end
