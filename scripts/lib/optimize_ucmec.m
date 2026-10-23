function out = optimize_ucmec(sc, scheme)
scheme = canonical_scheme(scheme);
switch scheme
    case 'proposed'
        out = admm_paper(sc);
    case {'oo', 'ro'}
        out = solve_clustered(sc, scheme);
    case 'bcd'
        out = solve_bcd(sc);
    case 'so'
        out = solve_so(sc);
    otherwise
        error('Unknown scheme "%s".', scheme);
end
out.scheme = scheme;
end

function scheme = canonical_scheme(scheme)
scheme = lower(strtrim(scheme));
switch scheme
    case {'proposed'}
        scheme = 'proposed';
    case {'oo', 'optimal', 'optimal_offloading'}
        scheme = 'oo';
    case {'ro', 'raft'}
        scheme = 'ro';
    case {'so', 'single', 'single_ap'}
        scheme = 'so';
    case {'bcd', 'bcdo'}
        scheme = 'bcd';
    otherwise
        error('Unknown scheme "%s".', scheme);
end
end

function out = solve_clustered(sc, scheme)
% Paper Alg.2 ADMM: continuous a in [0,1]; feasible clustering is local copy Â
% on the simplex (sum_m â_{m,n}=1). Stop when ||A-Â||_F <= gamma.
[A_hat, ~, ~] = initial_point(sc);
A = A_hat;
lambda = zeros(sc.M, sc.N);
if strcmp(scheme, 'ro')
    P_fixed = rand(sc.M, 1);
else
    P_fixed = [];
end
reserve = consensus_reserve(sc);
verbose = ~isfield(sc, 'verbose') || sc.verbose;
q_eff = sc.q * 4e-4;
use_admm = strcmp(scheme, 'proposed');
if use_admm
    tmax = max(sc.tmax, 100);
    a_iters = 40;
elseif strcmp(scheme, 'oo')
    % Weaker search than Proposed (paper: OO is a baseline, not full ADMM).
    tmax = 8;
    a_iters = 12;
    reserve = min(reserve, 0.1 * sc.B_MHz(:));
else
    tmax = 12;
    a_iters = 20;
end

[metrics, B, C] = assess(sc, A_hat, scheme, P_fixed, reserve);
history = empty_history(tmax + 1);
t = 1;
history = push_history(history, t, metrics, inf);
bestA = A_hat;
bestB = B;
bestC = C;
bestMet = metrics;
bestValue = objective_of(metrics, scheme);
if verbose
    fprintf('  [%s] init     Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
        scheme, metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
end

for sweep = 1:tmax
    % Global update (P4/P5 style): continuous a_{m,n} in [0,1], ADMM-coupled to Â.
    order = user_order(metrics);
    for n = order
        A = improve_user_continuous(sc, A, n, A_hat, lambda, q_eff, scheme, P_fixed, reserve, a_iters);
    end
    if use_admm
        % Local update (P6'): Â on simplex.
        for n = 1:sc.N
            A_hat(:, n) = project_simplex(A(:, n) + lambda(:, n) / max(q_eff, 1e-12));
        end
        gap = norm(A - A_hat, 'fro');
        lambda = lambda + q_eff * (A - A_hat);
        A_eval = A_hat;
    else
        % OO/RO: continuous a, then project to feasible clustering.
        for n = 1:sc.N
            A(:, n) = project_simplex(max(A(:, n), 0));
        end
        gap = norm(A - A_hat, 'fro');
        A_hat = A;
        A_eval = A;
    end
    [metrics, B, C] = assess(sc, A_eval, scheme, P_fixed, reserve);
    t = t + 1;
    history = push_history(history, t, metrics, gap);
    if verbose
        fprintf('  [%s] sweep %2d  gap %.3e  Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
            scheme, sweep, gap, metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
    end
    if objective_of(metrics, scheme) < bestValue
        bestValue = objective_of(metrics, scheme);
        bestA = A_eval;
        bestB = B;
        bestC = C;
        bestMet = metrics;
    end
    if gap <= sc.gamma
        break
    end
end

[bestMet, bestB, bestC] = finish_cvx(sc, bestA, scheme, P_fixed);
out.A = bestA;
out.B = bestB;
out.C = bestC;
out.metrics = bestMet;
out.history = trim_history(history, t);
out.iterations = t;
end

function out = solve_bcd(sc)
% Centralized BCD (paper BCDO): fewer outer passes than ADMM → weaker at large N.
A = initial_association(sc);
reserve = consensus_reserve(sc);
verbose = ~isfield(sc, 'verbose') || sc.verbose;
n_outer = 3;
a_iters = 12;
history = empty_history(n_outer + 1);
[metrics, B, C] = assess(sc, A, 'proposed', [], reserve);
history = push_history(history, 1, metrics, inf);
lambda0 = zeros(size(A));
for pass = 1:n_outer
    A_prev = A;
    for n = user_order(metrics)
        A = improve_user_continuous(sc, A, n, A, lambda0, 0, 'proposed', [], reserve, a_iters);
        A(:, n) = project_simplex(max(A(:, n), 0));
    end
    [metrics, B, C] = assess(sc, A, 'proposed', [], reserve);
    gap = norm(A - A_prev, 'fro');
    history = push_history(history, pass + 1, metrics, gap);
    if gap <= sc.gamma
        break
    end
end
if verbose
    fprintf('  [bcd] Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
        metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
end
[metrics, B, C] = finish_cvx(sc, A, 'proposed', []);
out.A = A;
out.B = B;
out.C = C;
out.metrics = metrics;
out.history = trim_history(history, pass + 1);
out.iterations = pass + 1;
end

function A = initial_association(sc)
A = zeros(sc.M, sc.N);
for n = 1:sc.N
    [~, order] = sort(sc.beta_user(:, n), 'descend');
    A(order(1), n) = 1;
end
end

function out = solve_so(sc)
A = initial_association(sc);
reserve = consensus_reserve(sc);
verbose = ~isfield(sc, 'verbose') || sc.verbose;
[metrics, B, C] = assess(sc, A, 'proposed', [], reserve);
if verbose
    fprintf('  [so] Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
        metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
end
[metrics, B, C] = finish_cvx(sc, A, 'so', []);
out.A = A;
out.B = B;
out.C = C;
out.metrics = metrics;
out.history = [];
out.iterations = 1;
end

function order = user_order(metrics)
[~, order] = sort(metrics.D_user, 'descend');
end

function A = improve_user_continuous(sc, A, n, anchor, lambda, q_eff, scheme, P_fixed, reserve, a_iters)
% Continuous a_{:,n} in [0,1]^M (paper); simplex is enforced on local Â.
if nargin < 10 || isempty(a_iters)
    a_iters = 40;
end
M = sc.M;
a0 = min(max(A(:, n), 0), 1);
if sum(a0) <= 1e-12
    a0 = ones(M, 1) / M;
end
obj = @(a) continuous_score(sc, A, n, a, anchor, lambda, q_eff, scheme, P_fixed, reserve);
opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxIterations', a_iters, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8);
try
    a = fmincon(obj, a0, [], [], [], [], zeros(M, 1), ones(M, 1), [], opts);
catch
    a = projected_gradient_box(obj, a0, max(8, round(a_iters / 2)));
end
A(:, n) = min(max(a(:), 0), 1);
end

function score = continuous_score(sc, A, n, a, anchor, lambda, q_eff, scheme, P_fixed, reserve)
Trial = A;
Trial(:, n) = min(max(a(:), 0), 1);
% Physical rates/energy use a feasible (simplex) clustering for this user.
Trial_eval = Trial;
Trial_eval(:, n) = project_simplex(Trial(:, n));
met = assess(sc, Trial_eval, scheme, P_fixed, reserve);
score = objective_of(met, scheme) + admm_penalty(Trial, anchor, lambda, q_eff) ...
    + 40 * sum(max(0, met.D_user - sc.Dt));
end

function a = projected_gradient_box(obj, a0, n_iter)
a = min(max(a0(:), 0), 1);
f0 = obj(a);
eps_fd = 1e-4;
step = 0.25;
for it = 1:n_iter
    g = zeros(size(a));
    for i = 1:numel(a)
        e = zeros(size(a));
        e(i) = eps_fd;
        g(i) = (obj(min(max(a + e, 0), 1)) - f0) / eps_fd;
    end
    a_try = min(max(a - step * g, 0), 1);
    f_try = obj(a_try);
    if f_try < f0
        a = a_try;
        f0 = f_try;
        step = min(step * 1.2, 1);
    else
        step = step * 0.5;
        if step < 1e-5
            break
        end
    end
end
end

function value = objective_of(met, scheme)
if strcmp(scheme, 'oo')
    value = met.E_off;
else
    value = met.E_tot;
end
end

function p = admm_penalty(A, anchor, lambda, q_eff)
D = A - anchor;
p = 0.5 * q_eff * sum(D(:).^2) + sum(lambda(:) .* D(:));
end

function a = project_simplex(v)
% Euclidean projection onto {a >= 0, sum(a)=1}.
v = v(:);
n = numel(v);
if n == 0
    a = v;
    return
end
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

function [metrics, B, C] = assess(sc, A, scheme, P_fixed, reserve)
se_user = uplink_se(sc, A);
if use_cvx_in_search(sc)
    [metrics, B, C] = allocate_cvx(sc, A, se_user, scheme, P_fixed);
else
    [B, C] = allocate_resources(sc, A, se_user, scheme, P_fixed, reserve);
    metrics = evaluate_solution(sc, A, B, C, se_user, scheme, P_fixed);
end
end

function tf = use_cvx_in_search(sc)
tf = (isfield(sc, 'cvx_in_search') && sc.cvx_in_search) ...
    || (isappdata(0, 'ucmec_cvx_in_search') && getappdata(0, 'ucmec_cvx_in_search'));
end

function [metrics, B, C] = finish_cvx(sc, A, scheme, P_fixed)
se_user = uplink_se(sc, A);
[metrics, B, C] = allocate_cvx(sc, A, se_user, scheme, P_fixed);
end

function [B, C] = allocate_resources(sc, A, se_user, scheme, P_fixed, reserve)
M = sc.M;
N = sc.N;
use_consensus = ~strcmp(scheme, 'oo');
is_leader = false(M, 1);
if strcmp(scheme, 'ro')
    is_leader(elected_leader(sc, sc.B_MHz(:), sc.C_GHz(:), P_fixed)) = true;
elseif use_consensus
    is_leader(1) = true;
end
if ~use_consensus
    reserve = min(reserve(:), 0.1 * sc.B_MHz(:));
end
B = zeros(M, N);
C = zeros(M, N);
for pass = 1:4
    for m = 1:M
        active = A(m, :) > 1e-8;
        [B(m, :), C(m, :)] = waterfill_ap(sc, m, A(m, :), se_user, is_leader(m), ...
            reserve(m), use_consensus, active);
    end
    if ~use_consensus
        break
    end
    if strcmp(scheme, 'ro')
        continue
    end
    [B_bar, C_bar] = residuals(sc, B, C);
    is_leader = false(M, 1);
    is_leader(elected_leader(sc, B_bar, C_bar, [])) = true;
end
end

function [b, c] = waterfill_ap(sc, m, a_row, se_user, is_leader, reserve_m, use_consensus, active)
N = sc.N;
b = zeros(1, N);
c = zeros(1, N);
if ~any(active)
    return
end
a = a_row(:);
se = max(se_user(:), 1e-3);
idx = find(active);
Btot = sc.B_MHz(m);
Ctot = sc.C_GHz(m);
reserve_m = min(max(reserve_m, 1e-4), 0.85 * Btot);
C_floor = min(0.4 * Ctot, max((sc.Lb / 1e6) / 0.005, 1e-3));
C_budget = max(Ctot - C_floor, 0.5 * Ctot);

K = a .* sc.L(:) .* sc.rho(:) / 1e6;
U = sc.L(:) ./ (se * 1e3);
U(a <= 1e-8) = 0;
u = sc.eps_c * a .* U;
cp_share = max(sc.Dt(:), 1e-3);
bmin = zeros(N, 1);
cmin = zeros(N, 1);
bmin(idx) = a(idx) .* U(idx) ./ cp_share(idx);
cmin(idx) = K(idx) ./ cp_share(idx);

B_bar = max(reserve_m, 0.25 * Btot);
Kb = 0;
for it = 1:5
    if use_consensus
        [T_lead, T_follow, T_gen] = node_times(sc, m, B_bar, Ctot);
        if is_leader
            Kb = sc.eps_c * max(T_lead - T_gen, 0) * B_bar;
        else
            Kb = sc.eps_c * T_follow * B_bar;
        end
    end
    [b_vec, B_bar] = waterfill_bandwidth(u, Kb, Btot, reserve_m, bmin);
end
c_vec = waterfill_cpu(K, C_budget, cmin);
b = b_vec.';
c = c_vec.';
end

function [b, B_bar] = waterfill_bandwidth(u, Kb, Btot, reserve_m, bmin)
N = numel(u);
b = zeros(N, 1);
active = u > 0;
if ~any(active)
    B_bar = Btot;
    return
end
bmin = max(bmin(:), 0);
floor_sum = sum(bmin(active));
room = Btot - reserve_m;
if floor_sum >= room
    b(active) = bmin(active) * room / floor_sum;
    B_bar = reserve_m;
    return
end
extra_budget = Btot - floor_sum;
uu = u(active);
den = sum(sqrt(uu)) + sqrt(max(Kb, 0));
if den <= 0
    B_bar = reserve_m;
    b(active) = bmin(active) + (room - floor_sum) / nnz(active);
    return
end
B_bar = sqrt(max(Kb, 0)) * extra_budget / den;
shares = sqrt(uu) * extra_budget / den;
if B_bar < reserve_m || Kb <= 0
    B_bar = reserve_m;
    rest = room - floor_sum;
    shares = sqrt(uu) / max(sum(sqrt(uu)), 1e-12) * rest;
end
b(active) = bmin(active) + shares;
s = sum(b);
if s + B_bar > Btot
    overflow = s + B_bar - Btot;
    b(active) = bmin(active) + max(0, shares - overflow * shares / max(sum(shares), 1e-12));
    B_bar = Btot - sum(b);
end
B_bar = max(B_bar, reserve_m);
if sum(b) + B_bar > Btot
    b(active) = b(active) * (Btot - B_bar) / max(sum(b), 1e-12);
end
end

function c = waterfill_cpu(K, C_budget, cmin)
c = zeros(size(K));
active = K > 0;
if ~any(active)
    return
end
cmin = max(cmin(:), 0);
if sum(cmin(active)) >= C_budget
    c(active) = cmin(active) * C_budget / sum(cmin(active));
    return
end
rest = C_budget - sum(cmin(active));
w = sqrt(K(active));
c(active) = cmin(active) + w / sum(w) * rest;
end

function metrics = evaluate_solution(sc, A, B, C, se_user, scheme, P_fixed)
[B_bar, C_bar] = residuals(sc, B, C);
if strcmp(scheme, 'ro')
    leader = elected_leader(sc, B_bar, C_bar, P_fixed);
else
    leader = elected_leader(sc, B_bar, C_bar, []);
end
se_user = max(se_user(:), 1e-3);
De = zeros(sc.M, sc.N);
tx_time = zeros(1, sc.N);
rate_user = zeros(1, sc.N);
for n = 1:sc.N
    band = 0;
    for m = 1:sc.M
        if A(m, n) <= 1e-8 || C(m, n) <= 0
            continue
        end
        De(m, n) = A(m, n) * sc.L(n) * sc.rho(n) / (C(m, n) * 1e6);
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
T_node = zeros(sc.M, 1);
for m = 1:sc.M
    [T_lead, T_follow, T_gen] = node_times(sc, m, B_bar(m), C_bar(m));
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

function [T_lead, T_follow, T_gen] = node_times(sc, m, B_bar, C_bar)
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
    mp = others(k);
    se_s(k) = link_se(sc, m, mp, Bs_pkt);
    se_b(k) = link_se(sc, m, mp, Bb);
end
n_link = numel(others);
T_state = n_link * rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * Bs * 1e6);
T_block = n_link * (sc.Lb * 1e3) / (mean(se_b) * Bb * 1e6);
T_lead = T_state + T_block + T_gen;
T_follow = rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * Bs * 1e6);
end

function se = link_se(sc, m, mp, B_mhz)
Bh = max(B_mhz, 1e-9) * 1e6;
sinr = sc.S_ap(m, mp) / (sc.N0 * Bh);
se = log2(1 + max(sinr, 0));
se = max(se, 1e-3);
end

function reserve = consensus_reserve(sc)
reserve = zeros(sc.M, 1);
for m = 1:sc.M
    B = 1;
    Cref = max(sc.C_GHz(m), 1e-3);
    for k = 1:8 %#ok<NASGU>
        [T_lead, ~, T_gen] = node_times(sc, m, B, Cref);
        tx = max(T_lead - T_gen, 1e-9);
        B = B * tx / max(sc.block_interval, 1e-6);
    end
    reserve(m) = min(max(B, 1e-4), 0.85 * sc.B_MHz(m));
end
end

function idx = elected_leader(sc, B_bar, C_bar, P_fixed)
if nargin >= 4 && ~isempty(P_fixed)
    [~, idx] = max(P_fixed(:));
    return
end
R = C_bar(:) / sc.Lb + max(B_bar(:), 0) / (sc.Ls + sc.Lb);
[~, idx] = max(R);
end

function [A, B, C] = initial_point(sc, A_fixed)
if nargin < 2 || isempty(A_fixed)
    A = initial_association(sc);
else
    A = A_fixed;
end
B = zeros(sc.M, sc.N);
C = zeros(sc.M, sc.N);
end

function [B_bar, C_bar] = residuals(sc, B, C)
B_bar = max(sc.B_MHz(:) - sum(B, 2), 1e-6);
C_bar = max(sc.C_GHz(:) - sum(C, 2), 1e-8);
end

function h = empty_history(n)
h.E_tot = nan(1, n);
h.E_off = nan(1, n);
h.E_con = nan(1, n);
h.D_off_ms = nan(1, n);
h.D_con_ms = nan(1, n);
h.gap = nan(1, n);
end

function h = push_history(h, t, metrics, gap)
h.E_tot(t) = metrics.E_tot;
h.E_off(t) = metrics.E_off;
h.E_con(t) = metrics.E_con;
h.D_off_ms(t) = metrics.D_off_ms;
h.D_con_ms(t) = metrics.D_con_ms;
h.gap(t) = gap;
end

function h = trim_history(h, t)
fields = fieldnames(h);
for k = 1:numel(fields)
    h.(fields{k}) = h.(fields{k})(1:t);
end
end
