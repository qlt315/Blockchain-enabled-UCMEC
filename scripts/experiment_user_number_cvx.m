function experiment_user_number_cvx()
% User-number sweep only, with CVX inside association search.
setappdata(0, 'ucmec_cvx_in_search', true);
cleaner = onCleanup(@() setappdata(0, 'ucmec_cvx_in_search', false)); %#ok<NASGU>
fprintf('CVX-in-search enabled for user_number sweep.\n');
experiment_user_number();
end
