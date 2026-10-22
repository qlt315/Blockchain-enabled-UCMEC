function energy = sweep_values(schemes, values, override_fn, seeds, label)
energy = zeros(numel(schemes), numel(values));
for s = 1:numel(seeds)
    for j = 1:numel(values)
        if s == 1
            fprintf('%s %g\n', label, values(j));
        end
        sc = build_scenario(5, 5, seeds(s), override_fn(values(j)));
        for i = 1:numel(schemes)
            energy(i, j) = energy(i, j) + optimize_ucmec(sc, schemes{i}).metrics.E_tot;
        end
    end
end
energy = energy / numel(seeds);
end
