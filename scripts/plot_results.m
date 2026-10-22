function plot_results()
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
repo = setup_path();
set(0, 'DefaultFigureVisible', 'off');

out_dir = fullfile(repo, 'results');
fig_dir = fullfile(out_dir, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

styles = {'-o', '-d', '-*', '-x', '-^'};
labels = {'OO', 'SO', 'RO', 'Proposed', 'BCDO'};

ap_file = fullfile(out_dir, 'convergence_ap.mat');
if exist(ap_file, 'file')
    S = load(ap_file);
    plot_convergence(S.convergence_ap, S.M_list, 30, styles, ...
        @(m) sprintf('Proposed, M=%d', m), fullfile(fig_dir, 'convergence_ap.png'));
end

q_file = fullfile(out_dir, 'convergence_q.mat');
if exist(q_file, 'file')
    S = load(q_file);
    plot_convergence(S.convergence_q, S.q_list, 40, styles, ...
        @(q) sprintf('Proposed, q=%g', q), fullfile(fig_dir, 'convergence_q.png'));
end

user_file = fullfile(out_dir, 'user_number.mat');
if exist(user_file, 'file')
    S = load(user_file);
    plot_group(S.N_list, S.rate, styles, labels, 'Number of users', ...
        'Average uplink rate (Mbps)', fullfile(fig_dir, 'rate_users.png'));
    plot_group(S.N_list, S.D_off, styles, labels, 'Number of users', ...
        'Average offloading delay (ms)', fullfile(fig_dir, 'delay_offloading.png'));
    plot_group(S.N_list, S.D_con, styles, labels, 'Number of users', ...
        'Consensus delay (ms)', fullfile(fig_dir, 'delay_consensus.png'));
    plot_group(S.N_list, S.D_tot, styles, labels, 'Number of users', ...
        'Total delay (ms)', fullfile(fig_dir, 'delay_total.png'));
    plot_group(S.N_list, S.E_off, styles, labels, 'Number of users', ...
        'Offloading energy consumption (J)', fullfile(fig_dir, 'energy_offloading.png'));
    plot_group(S.N_list, S.E_con, styles, labels, 'Number of users', ...
        'Consensus energy consumption (J)', fullfile(fig_dir, 'energy_consensus.png'));
    plot_group(S.N_list, S.E_tot, styles, labels, 'Number of users', ...
        'Total energy consumption (J)', fullfile(fig_dir, 'energy_total.png'));
end

plot_saved(out_dir, fig_dir, 'bandwidth.mat', 'B_list', 'E_B', styles, labels, ...
    'Bandwidth of APs (MHz)', 'Total energy consumption (J)', 'impact_bandwidth.png');
plot_saved(out_dir, fig_dir, 'cpu.mat', 'C_list', 'E_C', styles, labels, ...
    'Computing resources of APs (GHz)', 'Total energy consumption (J)', 'impact_cpu.png');
plot_saved(out_dir, fig_dir, 'epsilon.mat', 'ratio_list', 'E_ratio', styles, labels, ...
    'Ratio of \epsilon_c to \epsilon_p', 'Total energy consumption (J)', 'impact_epsilon.png');
plot_saved(out_dir, fig_dir, 'density.mat', 'rho_list', 'E_rho', styles, labels, ...
    'Computing density (CPU cycles/bit)', 'Total energy consumption (J)', 'impact_density.png');
plot_saved(out_dir, fig_dir, 'block_size.mat', 'Lb_list', 'E_Lb', styles, labels, ...
    'Size of block (kbits)', 'Total energy consumption (J)', 'impact_blocksize.png');
plot_saved(out_dir, fig_dir, 'ap_number.mat', 'M_list', 'E_M', styles, labels, ...
    'Number of APs', 'Total energy consumption (J)', 'impact_aps.png');
plot_saved(out_dir, fig_dir, 'ap_number.mat', 'M_list', 'rate_M', styles, labels, ...
    'Number of APs', 'Average uplink rate (Mbps)', 'rate_aps.png');
end

function plot_convergence(Y, keys, n_show, styles, label_fn, path)
figure('Color', 'w');
hold on;
n_show = min(n_show, size(Y, 2));
x = 0:n_show - 1;
for i = 1:size(Y, 1)
    y = Y(i, 1:n_show);
    y = cummin(y);
    plot(x, y, styles{i}, 'LineWidth', 1.6, 'MarkerSize', 5);
end
xlabel('Iteration index');
ylabel('Total energy consumption (J)');
legend(arrayfun(label_fn, keys, 'UniformOutput', false), 'Location', 'best');
style_axes();
saveas(gcf, path);
end

function plot_saved(out_dir, fig_dir, file_name, x_name, y_name, styles, labels, x_label, y_label, png_name)
file_path = fullfile(out_dir, file_name);
if ~exist(file_path, 'file')
    return
end
S = load(file_path);
plot_group(S.(x_name), S.(y_name), styles, labels, x_label, y_label, fullfile(fig_dir, png_name));
end

function plot_group(x, Y, styles, labels, x_label, y_label, path)
figure('Color', 'w');
hold on;
for i = 1:size(Y, 1)
    plot(x, Y(i, :), styles{i}, 'LineWidth', 2.2, 'MarkerSize', 8);
end
xlabel(x_label);
ylabel(y_label);
legend(labels, 'Location', 'best');
style_axes();
saveas(gcf, path);
end

function style_axes()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end
