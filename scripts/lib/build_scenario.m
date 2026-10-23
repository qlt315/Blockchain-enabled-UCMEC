function sc = build_scenario(M, N, seed, override)
if nargin < 3 || isempty(seed)
    seed = 1;
end
if nargin < 4 || isempty(override)
    override = struct();
end

Mpool = max(M, 12);
Npool = max(N, 14);
rng(seed);

sc.M = M;
sc.N = N;
sc.varrho = 2;
sc.radius = 200;
sc.P_user = 0.1;
sc.P_ap = 0.2;
sc.N0 = 10^((-174 - 30) / 10);
sc.eps_c = 1.5;
sc.eps_p = 1;
sc.Ls = 10;
sc.Lb = 500;
sc.Di = 1;
sc.block_interval = 1;
sc.q = 10;
sc.gamma = 0.01;
sc.tmax = 100;
sc.verbose = false;
sc.cvx_in_search = false;
sc.seed = seed;

names = fieldnames(override);
for k = 1:numel(names)
    sc.(names{k}) = override.(names{k});
end

theta = 2 * pi * rand(Mpool, 1);
rad = sc.radius * sqrt(rand(Mpool, 1));
x_ap = rad .* cos(theta);
y_ap = rad .* sin(theta);
theta = 2 * pi * rand(Npool, 1);
rad = sc.radius * sqrt(rand(Npool, 1));
x_user = rad .* cos(theta);
y_user = rad .* sin(theta);
sc.x_ap = x_ap(1:M);
sc.y_ap = y_ap(1:M);
sc.x_user = x_user(1:N);
sc.y_user = y_user(1:N);

B_pool = 10 + 40 * rand(1, Mpool);
C_pool = 5 + 15 * rand(1, Mpool);
L_pool = 100 + 100 * rand(1, Npool);
rho_pool = 1000 + 1000 * rand(1, Npool);
Dt_pool = (100 + 50 * rand(1, Npool)) * 1e-3;

if isfield(override, 'B_MHz')
    B_pool = override.B_MHz(:)';
    if isscalar(B_pool)
        B_pool = B_pool * ones(1, Mpool);
    end
end
if isfield(override, 'C_GHz')
    C_pool = override.C_GHz(:)';
    if isscalar(C_pool)
        C_pool = C_pool * ones(1, Mpool);
    end
end
if isfield(override, 'L_kbit')
    L_pool = override.L_kbit(:)';
    if isscalar(L_pool)
        L_pool = L_pool * ones(1, Npool);
    end
end
if isfield(override, 'rho')
    rho_pool = override.rho(:)';
    if isscalar(rho_pool)
        rho_pool = rho_pool * ones(1, Npool);
    end
end
if isfield(override, 'Dt_ms')
    Dt_pool = override.Dt_ms(:)' * 1e-3;
    if isscalar(Dt_pool)
        Dt_pool = Dt_pool * ones(1, Npool);
    end
end
sc.B_MHz = B_pool(1:M);
sc.C_GHz = C_pool(1:M);
sc.L = L_pool(1:N);
sc.rho = rho_pool(1:N);
sc.Dt = Dt_pool(1:N);

sc.verbose = false;
if ~isfield(override, 'tmax')
    sc.tmax = 100;
end
if ~isfield(override, 'q')
    sc.q = 10;
end

beta_user = zeros(Mpool, Npool);
for m = 1:Mpool
    for n = 1:Npool
        d_km = max(hypot(x_ap(m) - x_user(n), y_ap(m) - y_user(n)), 1) / 1e3;
        beta_user(m, n) = 10^(-pl_db(d_km) / 10);
    end
end
beta_ap = zeros(Mpool, Mpool);
ap_frob = zeros(Mpool, Mpool);
for m = 1:Mpool
    for mp = 1:Mpool
        if m == mp
            continue
        end
        d_km = max(hypot(x_ap(m) - x_ap(mp), y_ap(m) - y_ap(mp)), 1) / 1e3;
        beta_ap(m, mp) = 10^(-pl_db(d_km) / 10);
        H = (randn(sc.varrho, sc.varrho) + 1i * randn(sc.varrho, sc.varrho)) / sqrt(2);
        ap_frob(m, mp) = sum(abs(H(:)).^2);
    end
end
g = zeros(sc.varrho, Mpool, Npool);
for m = 1:Mpool
    for n = 1:Npool
        h = (randn(sc.varrho, 1) + 1i * randn(sc.varrho, 1)) / sqrt(2);
        g(:, m, n) = sqrt(beta_user(m, n)) * h;
    end
end

sc.beta_user = beta_user(1:M, 1:N);
sc.beta_ap = beta_ap(1:M, 1:M);
sc.ap_frob = ap_frob(1:M, 1:M);
sc.g = g(:, 1:M, 1:N);
sc.S_ap = backhaul_power(sc);
end

function pl = pl_db(d_km)
pl = 128.1 + 37.6 * log10(d_km);
end

function S = backhaul_power(sc)
M = sc.M;
S = zeros(M, M);
for m = 1:M
    for mp = 1:M
        if m == mp
            continue
        end
        S(m, mp) = sc.P_ap * sc.beta_ap(m, mp) * sc.ap_frob(m, mp);
    end
end
end
