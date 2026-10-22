function out = optimize_ucmec(sc, scheme)
scheme = canonical_scheme(scheme);
switch scheme
    case {'proposed', 'oo', 'ro'}
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
[A, ~, ~] = initial_point(sc);
lambda = zeros(sc.M, sc.N);
if strcmp(scheme, 'ro')
    P_fixed = rand(sc.M, 1);
else
    P_fixed = [];
end
reserve = consensus_reserve(sc);
verbose = ~isfield(sc, 'verbose') || sc.verbose;
q_eff = sc.q * 4e-4;
tmax = sc.tmax;
if strcmp(scheme, 'proposed')
    tmax = max(sc.tmax, 20);
    q_eff = sc.q * 2e-4;
elseif strcmp(scheme, 'oo')
    reserve = min(reserve, 0.08 * sc.B_MHz(:));
end

[metrics, B, C] = assess(sc, A, scheme, P_fixed, reserve);
history = empty_history(tmax + 1);
t = 1;
history = push_history(history, t, metrics, inf);
bestA = A;
bestB = B;
bestC = C;
bestMet = metrics;
bestValue = objective_of(metrics, scheme);
if verbose
    fprintf('  [%s] init     Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
        scheme, metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
end

for sweep = 1:tmax
    anchor = A;
    order = user_order(metrics);
    if strcmp(scheme, 'proposed')
        mode = 'rich';
    else
        mode = 'pair';
    end
    for n = order
        A = improve_user(sc, A, n, anchor, lambda, q_eff, scheme, P_fixed, reserve, mode);
    end
    gap = norm(A - anchor, 'fro');
    lambda = lambda + q_eff * (A - anchor);
    [metrics, B, C] = assess(sc, A, scheme, P_fixed, reserve);
    t = t + 1;
    history = push_history(history, t, metrics, gap);
    if verbose
        fprintf('  [%s] sweep %2d  gap %.3e  Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
            scheme, sweep, gap, metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
    end
    if objective_of(metrics, scheme) < bestValue
        bestValue = objective_of(metrics, scheme);
        bestA = A;
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
A = initial_association(sc);
reserve = consensus_reserve(sc);
verbose = ~isfield(sc, 'verbose') || sc.verbose;
[metrics, B, C] = assess(sc, A, 'proposed', [], reserve);
history = empty_history(2);
history = push_history(history, 1, metrics, inf);
anchor = A;
for n = user_order(metrics)
    A_try = improve_user(sc, A, n, anchor, zeros(size(A)), 0, 'proposed', [], reserve, 'pair');
    if norm(A_try(:, n) - A(:, n)) > 1e-8
        A = A_try;
        [metrics, B, C] = assess(sc, A, 'proposed', [], reserve);
    end
end
gap = norm(A - anchor, 'fro');
history = push_history(history, 2, metrics, gap);
if verbose
    fprintf('  [bcd] Eoff %.3f  Econ %.3f  Etot %.3f J  Doff %.2f ms  Dcon %.2f ms\n', ...
        metrics.E_off, metrics.E_con, metrics.E_tot, metrics.D_off_ms, metrics.D_con_ms);
end
[metrics, B, C] = finish_cvx(sc, A, 'proposed', []);
out.A = A;
out.B = B;
out.C = C;
out.metrics = metrics;
out.history = history;
out.iterations = 2;
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

function A = improve_user(sc, A, n, anchor, lambda, q_eff, scheme, P_fixed, reserve, mode)
cols = candidate_columns(sc, n, mode);
best_col = A(:, n);
best_score = inf;
for k = 1:numel(cols)
    Trial = A;
    Trial(:, n) = cols{k};
    met = assess(sc, Trial, scheme, P_fixed, reserve);
    score = objective_of(met, scheme) + admm_penalty(Trial, anchor, lambda, q_eff) ...
        + 40 * sum(max(0, met.D_user - sc.Dt));
    if score < best_score
        best_score = score;
        best_col = cols{k};
    end
end
A(:, n) = best_col;
end

function cols = candidate_columns(sc, n, mode)
cols = cell(1, sc.M);
for m = 1:sc.M
    a = zeros(sc.M, 1);
    a(m) = 1;
    cols{m} = a;
end
if strcmp(mode, 'single') || sc.M < 2
    return
end
[~, order] = sort(sc.beta_user(:, n), 'descend');
if strcmp(mode, 'rich')
    top = order(1:min(4, sc.M));
    fracs = [0.5 0.5; 0.7 0.3; 0.3 0.7; 0.9 0.1; 0.1 0.9; 0.6 0.4; 0.4 0.6];
else
    top = order(1:min(3, sc.M));
    fracs = [0.5 0.5; 0.7 0.3; 0.3 0.7; 0.9 0.1];
end
for i = 1:numel(top)
    for j = i + 1:numel(top)
        for f = 1:size(fracs, 1)
            a = zeros(sc.M, 1);
            a(top(i)) = fracs(f, 1);
            a(top(j)) = fracs(f, 2);
            cols{end + 1} = a; %#ok<AGROW>
        end
    end
end
if strcmp(mode, 'rich') && numel(top) >= 3
    trips = [1/3 1/3 1/3; 0.5 0.25 0.25; 0.6 0.2 0.2];
    for t = 1:size(trips, 1)
        a = zeros(sc.M, 1);
        a(top(1)) = trips(t, 1);
        a(top(2)) = trips(t, 2);
        a(top(3)) = trips(t, 3);
        cols{end + 1} = a; %#ok<AGROW>
    end
end
end

function value = objective_of(met, scheme)
if strcmp(scheme, 'oo')
    value = met.E_off;
elseif strcmp(scheme, 'proposed')
    value = met.E_tot;
else
    value = met.E_tot;
end
end

function p = admm_penalty(A, anchor, lambda, q_eff)
D = A - anchor;
p = 0.5 * q_eff * sum(D(:).^2) + sum(lambda(:) .* D(:));
end

function [metrics, B, C] = assess(sc, A, scheme, P_fixed, reserve)
se_user = uplink_se(sc, A);
[B, C] = allocate_resources(sc, A, se_user, scheme, P_fixed, reserve);
metrics = evaluate_solution(sc, A, B, C, se_user, scheme, P_fixed);
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
    reserve = min(reserve(:), 0.08 * sc.B_MHz(:));
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

function se_user = uplink_se(sc, A)
N = sc.N;
sinr = zeros(1, N);
Phi_of = cell(1, N);
for n = 1:N
    Phi_of{n} = dominant_aps(A(:, n), sc.M);
end
for n = 1:N
    Phi = Phi_of{n};
    g_n = stack_g(sc, Phi, n);
    intra = [];
    for v = 1:N
        if v == n
            continue
        end
        if isequal(Phi_of{v}, Phi)
            intra(end + 1) = v; %#ok<AGROW>
        end
    end
    overlap = 0.08;
    if isfield(sc, 'overlap')
        overlap = sc.overlap;
    end
    noise = sc.N0 * mean(sc.B_MHz(Phi)) * 1e6;
    dof = numel(g_n);
    keep = min(numel(intra), max(dof - 1, 0));
    cancelled = [];
    if keep > 0
        strength = zeros(1, numel(intra));
        for k = 1:numel(intra)
            strength(k) = norm(stack_g(sc, Phi, intra(k)))^2;
        end
        [~, order] = sort(strength, 'descend');
        cancelled = intra(order(1:keep));
    end
    R = noise * eye(dof);
    for w = 1:N
        if w == n || any(cancelled == w)
            continue
        end
        g_w = stack_g(sc, Phi, w);
        R = R + overlap * sc.P_user * (g_w * g_w');
    end
    xi = R \ g_n;
    xi = xi / max(norm(xi), 1e-15);
    signal = sc.P_user * abs(xi' * g_n)^2;
    interf = 0;
    for w = 1:N
        if w == n || any(cancelled == w)
            continue
        end
        g_w = stack_g(sc, Phi, w);
        interf = interf + overlap * sc.P_user * abs(xi' * g_w)^2;
    end
    sinr(n) = signal / max(interf + noise, 1e-30);
end
sinr = min(max(sinr, 0), 1e8);
se_user = max(log2(1 + sinr), 1e-3);
end

function Phi = dominant_aps(a_col, K)
[vals, order] = sort(a_col, 'descend');
K = min(K, numel(order));
Phi = order(1:K);
Phi = Phi(vals(1:K) >= 0.02);
if isempty(Phi)
    Phi = order(1);
end
Phi = sort(Phi);
end

function g = stack_g(sc, Phi, user)
g = zeros(sc.varrho * numel(Phi), 1);
for k = 1:numel(Phi)
    idx = (k - 1) * sc.varrho + (1:sc.varrho);
    g(idx) = sc.g(:, Phi(k), user);
end
end
