function [repo, schemes, n_seed] = setup_path()
persistent ready repo_root
if isempty(repo_root)
    here = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(here));
end
if isempty(ready) || ~ready
    addpath(fileparts(mfilename('fullpath')));
    addpath(fullfile(repo_root, 'scripts'));
    if exist('cvx_begin', 'file') ~= 2
        addpath(fullfile(repo_root, 'cvx'));
        cvx_setup;
    end
    ready = true;
end
repo = repo_root;
out_dir = fullfile(repo, 'results');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
schemes = {'oo', 'so', 'ro', 'proposed', 'bcd'};
n_seed = 8;
end
