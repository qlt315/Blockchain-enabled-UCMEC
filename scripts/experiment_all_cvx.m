function experiment_all_cvx()
% Full-suite run with CVX inside association search (very slow).
% Usage: cd scripts; experiment_all_cvx; plot_results
setappdata(0, 'ucmec_cvx_in_search', true);
cleaner = onCleanup(@() setappdata(0, 'ucmec_cvx_in_search', false)); %#ok<NASGU>
fprintf('CVX-in-search enabled for this MATLAB session sweep.\n');
experiment_all();
end
