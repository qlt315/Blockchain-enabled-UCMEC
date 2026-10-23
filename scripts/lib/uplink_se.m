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
