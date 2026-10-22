function experiment_block_size()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
Lb_list = 100:100:1000;
E_Lb = sweep_values(schemes, Lb_list, @(v) struct('Lb', v), 1:n_seed, 'Block size');
save(fullfile(repo, 'results', 'block_size.mat'), 'schemes', 'Lb_list', 'E_Lb');
end
