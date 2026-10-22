# Energy-Efficient Blockchain-enabled UCMEC

This is the updated version of the MATLAB evaluation code for [IEEE TCCN 10460318](https://ieeexplore.ieee.org/abstract/document/10460318).

## Repository layout

```
Blockchain-enabled-UCMEC/
├── README.md
├── .gitignore
├── cvx/                 # CVX
├── scripts/             # experiment drivers and plotter
│   ├── lib/             # shared solvers and scenario builder
│   ├── experiment_*.m
│   ├── experiment_all.m
│   └── plot_results.m
└── results/             # .mat outputs and figures/
    └── figures/
```

### Top level

| Path | Role |
| --- | --- |
| `cvx/` | Bundled CVX First run calls `cvx_setup` via `setup_path`. |
| `scripts/` | One script per figure sweep, plus `plot_results`. |
| `results/` | Saved `.mat` arrays and PNG figures under `results/figures/`. |

### `scripts/`

| File | Role |
| --- | --- |
| `experiment_convergence.m` | Convergence vs iteration (AP count and penalty `q`). |
| `experiment_user_number.m` | Sweep user count `N = 5:14` (energy, delay, rate). |
| `experiment_bandwidth.m` | Sweep AP bandwidth `10:10:50` MHz. |
| `experiment_cpu.m` | Sweep AP CPU `10:10:50` GHz. |
| `experiment_epsilon.m` | Sweep `eps_c / eps_p` with `eps_c = 1.5`. |
| `experiment_density.m` | Sweep task density `500:500:5000` cycles/bit. |
| `experiment_block_size.m` | Sweep block size `100:100:1000` kbit. |
| `experiment_ap_number.m` | Sweep AP count `M = 3:12` (energy and uplink rate). |
| `experiment_all.m` | Runs every sweep above in order. |
| `plot_results.m` | Reads `results/*.mat` and writes `results/figures/*.png`. |

### `scripts/lib/`

| File | Role |
| --- | --- |
| `setup_path.m` | Adds paths, sets up CVX, returns repo root / scheme list / seed count (`n_seed = 8`). |
| `build_scenario.m` | Builds a nested-RNG scenario for given `M`, `N`, and seed. |
| `optimize_ucmec.m` | Association search (waterfill) and scheme entry points (`oo`, `so`, `ro`, `proposed`, `bcd`). |
| `allocate_cvx.m` | Final convex bandwidth/CPU allocation (CVX + SDPT3). |
| `sweep_values.m` | Shared helper for parameter sweeps over seeds and schemes. |

### Schemes

| Name | Alias | Meaning |
| --- | --- | --- |
| `proposed` | — | Joint association and allocation |
| `oo` | `optimal_offloading` | Offloading-first objective |
| `ro` | `raft` | Same search; leader sampled once uniformly |
| `so` | `single_ap` | Each user fixed to its strongest AP |
| `bcd` | — | Same objective as `proposed`, one coordinate pass |

## How to run

**Requirements:** MATLAB R2021b+ and the bundled CVX under `cvx/`.

From the `scripts/` folder:

```matlab
experiment_all      % run every sweep
plot_results        % regenerate figures from results/*.mat
```

Or run a single sweep:

```matlab
experiment_user_number
experiment_bandwidth
% ...
plot_results
```

Each sweep averages 8 nested seeds. Outputs land in `results/`; figures in `results/figures/`. The results might be a little bit different from the original paper since we have modified the realization and evaluation process. The system model, algorithms and the parameters are still align with the paper.
