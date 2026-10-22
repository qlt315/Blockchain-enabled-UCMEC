function experiment_density()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
[repo, schemes, n_seed] = setup_path();
rho_list = 500:500:5000;
E_rho = sweep_values(schemes, rho_list, @(v) struct('rho', v), 1:n_seed, 'Density');
save(fullfile(repo, 'results', 'density.mat'), 'schemes', 'rho_list', 'E_rho');
end
