function [metrics, B, C] = allocate_cvx(sc, A, se_user, scheme, P_fixed)
M = sc.M;
N = sc.N;
se_user = max(se_user(:)', 1e-3);
[k_lead, k_follow, k_gen] = consensus_coeff(sc);
if nargin >= 5 && strcmp(scheme, 'ro') && ~isempty(P_fixed)
    [~, leader] = max(P_fixed(:));
    [metrics, B, C] = solve_allocation(sc, A, se_user, scheme, k_lead, k_follow, k_gen, leader);
    return
end
reputation = sc.C_GHz(:) / sc.Lb + sc.B_MHz(:) / (sc.Ls + sc.Lb);
[~, leader] = max(reputation);
[metrics, B, C] = solve_allocation(sc, A, se_user, scheme, k_lead, k_follow, k_gen, leader);
reputation = metrics.C_bar(:) / sc.Lb + max(metrics.B_bar(:), 0) / (sc.Ls + sc.Lb);
[~, elected] = max(reputation);
if elected ~= leader
    [metrics, B, C] = solve_allocation(sc, A, se_user, scheme, k_lead, k_follow, k_gen, elected);
end
end

function [metrics, B, C] = solve_allocation(sc, A, se_user, scheme, k_lead, k_follow, k_gen, leader)
M = sc.M;
N = sc.N;
tx_coef = sc.L(:)' ./ (se_user * 1e3);
oo_mode = strcmp(scheme, 'oo');
cvx_solver sdpt3
cvx_precision low
cvx_begin quiet
    variables b(M, N) c(M, N) Bbar(M) Cbar(M) slack(N)
    expression tx(1, N)
    expression cpu_time(1, N)
    for n = 1:N
        tx(n) = tx_coef(n) * inv_pos(sum(b(:, n)));
        first = true;
        for m = 1:M
            if A(m, n) <= 1e-8
                continue
            end
            term = (A(m, n) * sc.L(n) * sc.rho(n) / 1e6) * inv_pos(c(m, n));
            if first
                cpu_time(n) = term;
                first = false;
            else
                cpu_time(n) = max(cpu_time(n), term);
            end
        end
    end
    E_off = sc.eps_c * sum(tx) + sc.eps_p * sum(cpu_time);
    E_con = 0;
    for m = 1:M
        if m == leader
            E_con = E_con + sc.eps_c * k_lead(m) * inv_pos(Bbar(m)) ...
                + sc.eps_p * k_gen * inv_pos(Cbar(m));
        else
            E_con = E_con + sc.eps_c * k_follow(m) * inv_pos(Bbar(m));
        end
    end
    if oo_mode
        % Offloading-first; light consensus weight keeps SO worst on totals.
        minimize(E_off + 0.25 * E_con + 40 * sum(slack))
    else
        minimize(E_off + E_con + 40 * sum(slack))
    end
    subject to
        sum(b, 2) + Bbar == sc.B_MHz(:)
        sum(c, 2) + Cbar == sc.C_GHz(:)
        b >= 0
        c >= 0
        slack >= 0
        Bbar >= 0.05
        Cbar >= 0.05
        tx + cpu_time <= sc.Dt(:)' + slack(:)'
        for m = 1:M
            for n = 1:N
                if A(m, n) <= 1e-8
                    b(m, n) == 0
                    c(m, n) == 0
                else
                    c(m, n) >= 1e-4
                    b(m, n) >= 1e-4
                end
            end
        end
cvx_end

if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
    metrics = failed_metrics(sc);
    B = zeros(M, N);
    C = zeros(M, N);
    return
end
metrics = metrics_from_solution(sc, A, se_user, b, c, Bbar, Cbar, leader, k_lead, k_follow, k_gen);
B = b;
C = c;
end

function [k_lead, k_follow, k_gen] = consensus_coeff(sc)
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
        se_s(k) = link_se(sc, m, others(k), B_ref * frac_s / rounds);
        se_b(k) = link_se(sc, m, others(k), B_ref * (1 - frac_s));
    end
    k_state = n_link * rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * frac_s * 1e6);
    k_block = n_link * (sc.Lb * 1e3) / (mean(se_b) * (1 - frac_s) * 1e6);
    k_lead(m) = k_state + k_block;
    k_follow(m) = rounds^2 * (sc.Ls * 1e3) / (mean(se_s) * frac_s * 1e6);
end
k_gen = (sc.Lb * 1e3) / 1e9;
end

function se = link_se(sc, m, mp, B_mhz)
Bh = max(B_mhz, 1e-9) * 1e6;
sinr = sc.S_ap(m, mp) / (sc.N0 * Bh);
se = log2(1 + max(sinr, 0));
se = max(se, 1e-3);
end

function metrics = metrics_from_solution(sc, A, se_user, b, c, Bbar, Cbar, leader, k_lead, k_follow, k_gen)
b = full(max(b, 0));
c = full(max(c, 0));
Bbar = full(max(Bbar(:), 1e-6));
Cbar = full(max(Cbar(:), 1e-6));
band = sum(b, 1);
rate = band .* se_user;
tx = sc.L(:)' ./ max(rate * 1e3, 1e-12);
De = zeros(sc.M, sc.N);
for m = 1:sc.M
    for n = 1:sc.N
        if A(m, n) > 1e-8 && c(m, n) > 0
            De(m, n) = A(m, n) * sc.L(n) * sc.rho(n) / (c(m, n) * 1e6);
        end
    end
end
cpu_time = max(De, [], 1);
E_off = sc.eps_c * sum(tx) + sc.eps_p * sum(cpu_time);
T = zeros(sc.M, 1);
E_con = 0;
for m = 1:sc.M
    gen = k_gen / Cbar(m);
    if m == leader
        T(m) = k_lead(m) / Bbar(m) + gen;
        E_con = E_con + sc.eps_c * (k_lead(m) / Bbar(m)) + sc.eps_p * gen;
    else
        T(m) = k_follow(m) / Bbar(m);
        E_con = E_con + sc.eps_c * T(m);
    end
end
D_user = full(tx + cpu_time);
metrics.E_off = full(E_off);
metrics.E_con = full(E_con);
metrics.E_tot = full(E_off + E_con);
metrics.D_user = D_user;
metrics.D_off_ms = full(mean(D_user) * 1e3);
metrics.D_con_ms = full(mean(T) * 1e3);
metrics.D_tot_ms = full(metrics.D_off_ms + metrics.D_con_ms);
metrics.rate_Mbps = full(mean(rate));
metrics.B_bar = Bbar;
metrics.C_bar = Cbar;
metrics.P_leader = leader;
metrics.mean_B_bar = mean(Bbar);
end

function metrics = failed_metrics(sc)
metrics.E_off = 1e3;
metrics.E_con = 1e3;
metrics.E_tot = 2e3;
metrics.D_user = 10 * ones(1, sc.N);
metrics.D_off_ms = 1e4;
metrics.D_con_ms = 1e4;
metrics.D_tot_ms = 2e4;
metrics.rate_Mbps = 0;
metrics.B_bar = 1e-3 * ones(sc.M, 1);
metrics.C_bar = 1e-3 * ones(sc.M, 1);
metrics.P_leader = 1;
metrics.mean_B_bar = 1e-3;
end
