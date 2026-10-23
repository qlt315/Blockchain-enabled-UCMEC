function out = admm_paper(sc)
% Paper Alg.2 ADMM with Alg.1 (CVX P5') global updates and CVX P6' local update.
M = sc.M;
N = sc.N;
q = sc.q;
gamma = sc.gamma;
tmax = sc.tmax;
verbose = ~isfield(sc, 'verbose') || sc.verbose;
zeta1 = 0.5;
zeta2 = 0.5;
n_dc = 4;
delta = 1e-2;
min_sweeps = min(10, tmax);

A_hat = initial_association_soft(sc);
lambda = zeros(M, N);
[k_lead, k_follow, k_gen] = consensus_coeff_local(sc);

A = A_hat;
B = zeros(M, N);
C = zeros(M, N);
history = empty_history_local(tmax + 1);
se = uplink_se(sc, A_hat);
[metrics, ~, ~] = allocate_cvx(sc, A_hat, se, 'proposed', []);
t = 1;
history = push_history_local(history, t, metrics, inf);
bestA = A_hat;
bestB = B;
bestC = C;
bestMet = metrics;
bestE = metrics.E_tot;

if verbose
    fprintf('  [proposed/alg2] init  Etot %.3f\n', metrics.E_tot);
end

for sweep = 1:tmax
    se = uplink_se(sc, A_hat);
    for m = 1:M
        [A(m, :), B(m, :), C(m, :)] = alg1_update_ap(sc, m, A_hat(m, :)', lambda(m, :)', ...
            q, se(:), zeta1, zeta2, n_dc, delta, k_follow(m), k_gen);
    end
    A_hat = p6_update_ahat(A, lambda, q);
    gap = norm(A - A_hat, 'fro');
    lambda = lambda + q * (A - A_hat);

    se = uplink_se(sc, A_hat);
    [metrics, Bf, Cf] = allocate_cvx(sc, A_hat, se, 'proposed', []);
    % Prefer Alg.1 bandwidth/CPU when feasible; otherwise CVX polish.
    if all(B(:) >= 0) && sum(B(:)) > 1e-6
        metrics_try = evaluate_from_BC(sc, A_hat, B, C, se);
        if metrics_try.E_tot < metrics.E_tot
            metrics = metrics_try;
        else
            B = Bf;
            C = Cf;
        end
    else
        B = Bf;
        C = Cf;
    end

    t = t + 1;
    history = push_history_local(history, t, metrics, gap);
    if verbose
        fprintf('  [proposed/alg2] sweep %2d  gap %.3e  Etot %.3f  Eoff %.3f  Econ %.3f\n', ...
            sweep, gap, metrics.E_tot, metrics.E_off, metrics.E_con);
    end
    if metrics.E_tot < bestE
        bestE = metrics.E_tot;
        bestA = A_hat;
        bestB = B;
        bestC = C;
        bestMet = metrics;
    end
    if sweep >= min_sweeps && gap <= gamma
        break
    end
end

[bestMet, bestB, bestC] = allocate_cvx(sc, bestA, uplink_se(sc, bestA), 'proposed', []);
out.A = bestA;
out.B = bestB;
out.C = bestC;
out.metrics = bestMet;
out.history = trim_history_local(history, t);
out.iterations = t;
out.scheme = 'proposed';
end

function [a_row, b_row, c_row] = alg1_update_ap(sc, m, ahat, lam, q, se, zeta1, zeta2, n_dc, delta, k_f, k_gen) %#ok<INUSD>
N = sc.N;
Bm = sc.B_MHz(m);
Cm = sc.C_GHz(m);
L = sc.L(:);
rho = sc.rho(:);
Dt = sc.Dt(:);
se = max(se(:), 1e-3);
ahat = ahat(:);
lam = lam(:);

b_floor = 0.05;
c_floor = 0.05;
b = max((Bm - b_floor) / max(N, 1), 1e-3) * ones(N, 1);
c = max((Cm - c_floor) / max(N, 1), 1e-3) * ones(N, 1);
bp = 1 ./ b;
cp = 1 ./ c;
a = min(max(ahat, 0), 1);
kappa = max(c_floor, 1e-3) / max(b_floor, 1e-3);

for it = 1:n_dc
    g_bp = -zeta2 * (Cm / Bm) ./ max(bp, 1e-6).^2;
    g_cp = -(zeta1 + zeta2) ./ max(cp, 1e-6).^2;
    g_k = zeta1 * Bm;

    cvx_solver sdpt3
    cvx_precision low
    cvx_begin quiet
        variables a_v(N) bp_v(N) cp_v(N) chi_v(N) psi_v(N) kappa_v
        E_off = 0;
        for n = 1:N
            % Paper (27a) uses 1/M share; keep it and rely on ADMM + trust region.
            E_off = E_off + (sc.eps_c / sc.M) * (L(n) / se(n)) * chi_v(n) ...
                + (sc.eps_p / sc.M) * (L(n) * rho(n)) * psi_v(n);
        end
        sumb = sum(inv_pos(bp_v));
        sumc = sum(inv_pos(cp_v));
        Bbar = Bm - sumb;
        Cbar = Cm - sumc;
        E_con = sc.eps_c * k_f * inv_pos(Bbar) + sc.eps_p * k_gen * inv_pos(Cbar);
        admm = lam' * (a_v - ahat) + (q / 2) * sum_square(a_v - ahat);
        minimize(E_off + E_con + admm - g_bp' * bp_v - g_cp' * cp_v - g_k * kappa_v)
        subject to
            0 <= a_v <= 1
            bp_v >= 1 / Bm
            cp_v >= 1 / Cm
            chi_v >= 0
            psi_v >= 0
            kappa_v >= 0
            sumb <= Bm - b_floor
            sumc <= Cm - c_floor
            Bbar >= b_floor
            Cbar >= c_floor
            % DC trust region (Alg.1); wide enough for multi-step ADMM progress
            sum_square(a_v - ahat) <= 2.0
            for n = 1:N
                (L(n) / se(n)) * chi_v(n) + L(n) * rho(n) * psi_v(n) <= Dt(n) + 0.05
                % McCormick / RLT (31a,b) (32a,b)
                chi_v(n) >= a_v(n) / Bm
                chi_v(n) <= bp_v(n) + a_v(n) / Bm - 1 / Bm
                psi_v(n) >= a_v(n) / Cm
                psi_v(n) <= cp_v(n) + a_v(n) / Cm - 1 / Cm
            end
            kappa_v <= Cm / b_floor
    cvx_end

    if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
        break
    end
    a = max(full(a_v), 0);
    bp = full(bp_v);
    cp = full(cp_v);
    kappa = full(kappa_v); %#ok<NASGU>
end

a_row = a(:)';
b_row = (1 ./ max(bp(:), 1e-6))';
c_row = (1 ./ max(cp(:), 1e-6))';
mask = a_row < 1e-6;
b_row(mask) = 0;
c_row(mask) = 0;
end

function A_hat = p6_update_ahat(A, lambda, q)
% Paper P6' via CVX (equivalent to simplex projection of A + lambda/q).
[M, N] = size(A);
target = A + lambda / max(q, 1e-12);
cvx_solver sdpt3
cvx_precision low
cvx_begin quiet
    variable A_hat(M, N)
    minimize(sum(sum_square(A_hat - target)))
    subject to
        A_hat >= 0
        A_hat <= 1
        sum(A_hat, 1) == 1
cvx_end
if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
    A_hat = zeros(M, N);
    for n = 1:N
        A_hat(:, n) = project_simplex_local(target(:, n));
    end
else
    A_hat = full(A_hat);
end
end

function A = initial_association_soft(sc)
% Feasible Â on simplex: top-2 APs (better warm start than single-AP).
A = zeros(sc.M, sc.N);
for n = 1:sc.N
    [~, order] = sort(sc.beta_user(:, n), 'descend');
    if sc.M >= 2
        A(order(1), n) = 0.7;
        A(order(2), n) = 0.3;
    else
        A(order(1), n) = 1;
    end
end
end

function a = project_simplex_local(v)
v = v(:);
n = numel(v);
u = sort(v, 'descend');
cssv = cumsum(u);
rho = find(u .* (1:n)' > (cssv - 1), 1, 'last');
if isempty(rho)
    a = ones(n, 1) / n;
    return
end
theta = (cssv(rho) - 1) / rho;
a = max(v - theta, 0);
s = sum(a);
if s <= 0
    a = ones(n, 1) / n;
else
    a = a / s;
end
end

function [k_lead, k_follow, k_gen] = consensus_coeff_local(sc)
M = sc.M;
rounds = 2 + sc.Di;
frac_s = sc.Ls / (sc.Ls + sc.Lb);
B_ref = 10;
k_lead = zeros(M, 1);
k_follow = zeros(M, 1);
n_link = max(M - 1, 1);
for m = 1:M
    others = [1:m - 1, m + 1:M];
    if isempty(others)
        continue
    end
    se_s = zeros(1, numel(others));
    se_b = zeros(1, numel(others));
    for k = 1:numel(others)
        se_s(k) = link_se_local(sc, m, others(k), B_ref * frac_s / rounds);
        se_b(k) = link_se_local(sc, m, others(k), B_ref * (1 - frac_s));
    end
    k_state = n_link * rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * frac_s * 1e6);
    k_block = n_link * (sc.Lb * 1e3) / (mean(se_b) * (1 - frac_s) * 1e6);
    k_lead(m) = k_state + k_block;
    k_follow(m) = rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * frac_s * 1e6);
end
k_gen = (sc.Lb * 1e3) / 1e9;
end

function se = link_se_local(sc, m, mp, B_mhz)
Bh = max(B_mhz, 1e-9) * 1e6;
sinr = sc.S_ap(m, mp) / (sc.N0 * Bh);
se = max(log2(1 + max(sinr, 0)), 1e-3);
end

function metrics = evaluate_from_BC(sc, A, B, C, se_user)
% Thin wrapper: reuse waterfill-path evaluator through temporary assess-like call.
reserve = max(sc.B_MHz(:) - sum(B, 2), 0.05);
% Build metrics consistently with allocate path units.
se_user = max(se_user(:)', 1e-3);
M = sc.M;
N = sc.N;
B_bar = max(sc.B_MHz(:) - sum(B, 2), 1e-6);
C_bar = max(sc.C_GHz(:) - sum(C, 2), 1e-8);
[~, ~, k_gen] = consensus_coeff_local(sc);
% leader by reputation
rep = C_bar / sc.Lb + B_bar / (sc.Ls + sc.Lb);
[~, leader] = max(rep);
De = zeros(M, N);
tx_time = zeros(1, N);
rate_user = zeros(1, N);
for n = 1:N
    band = 0;
    for m = 1:M
        if A(m, n) > 1e-8 && C(m, n) > 0
            De(m, n) = A(m, n) * sc.L(n) * sc.rho(n) / (C(m, n) * 1e6);
        end
        if B(m, n) > 0
            band = band + B(m, n);
        end
    end
    rate_user(n) = band * se_user(n);
    tx_time(n) = sc.L(n) / max(rate_user(n) * 1e3, 1e-12);
end
cpu_time = max(De, [], 1);
D_user = tx_time + cpu_time;
E_con = 0;
T_node = zeros(M, 1);
for m = 1:M
    [T_lead, T_follow, T_gen] = node_times_local(sc, m, B_bar(m), C_bar(m));
    if m == leader
        T_node(m) = T_lead;
        E_con = E_con + sc.eps_c * max(T_lead - T_gen, 0) + sc.eps_p * T_gen;
    else
        T_node(m) = T_follow;
        E_con = E_con + sc.eps_c * T_follow;
    end
end
metrics.E_off = sc.eps_c * sum(tx_time) + sc.eps_p * sum(cpu_time);
metrics.E_con = E_con;
metrics.E_tot = metrics.E_off + metrics.E_con;
metrics.D_user = D_user;
metrics.D_off_ms = mean(D_user) * 1e3;
metrics.D_con_ms = mean(T_node) * 1e3;
metrics.D_tot_ms = metrics.D_off_ms + metrics.D_con_ms;
metrics.rate_Mbps = mean(rate_user);
metrics.P_leader = leader;
metrics.mean_B_bar = mean(B_bar);
end

function [T_lead, T_follow, T_gen] = node_times_local(sc, m, B_bar, C_bar)
B_bar = max(B_bar, 1e-6);
C_bar = max(C_bar, 1e-8);
T_gen = (sc.Lb * 1e3) / (C_bar * 1e9);
if sc.M < 2
    T_lead = T_gen;
    T_follow = 0;
    return
end
frac_s = sc.Ls / (sc.Ls + sc.Lb);
Bs = max(B_bar * frac_s, 1e-6);
Bb = max(B_bar * (1 - frac_s), 1e-6);
rounds = 2 + sc.Di;
Bs_pkt = Bs / rounds;
others = [1:m - 1, m + 1:sc.M];
se_s = zeros(1, numel(others));
se_b = zeros(1, numel(others));
for k = 1:numel(others)
    se_s(k) = link_se_local(sc, m, others(k), Bs_pkt);
    se_b(k) = link_se_local(sc, m, others(k), Bb);
end
n_link = numel(others);
T_state = n_link * rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * Bs * 1e6);
T_block = n_link * (sc.Lb * 1e3) / (mean(se_b) * Bb * 1e6);
T_lead = T_state + T_block + T_gen;
T_follow = rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * Bs * 1e6);
end

function h = empty_history_local(cap)
h.E_tot = nan(1, cap);
h.E_off = nan(1, cap);
h.E_con = nan(1, cap);
h.D_off_ms = nan(1, cap);
h.D_con_ms = nan(1, cap);
h.gap = nan(1, cap);
end

function h = push_history_local(h, t, metrics, gap)
h.E_tot(t) = metrics.E_tot;
h.E_off(t) = metrics.E_off;
h.E_con(t) = metrics.E_con;
h.D_off_ms(t) = metrics.D_off_ms;
h.D_con_ms(t) = metrics.D_con_ms;
h.gap(t) = gap;
end

function h = trim_history_local(h, t)
fields = fieldnames(h);
for k = 1:numel(fields)
    h.(fields{k}) = h.(fields{k})(1:t);
end
end
