# Publication Plotting

For the full project run order, including when to generate plot data, see
`docs/GETTING_STARTED.md`.

`publication/plots` contains the legacy publication-figure live scripts. Those
scripts only load plot-data MAT files and make figures; they do not regenerate
data, load saved networks, or change plotting dependencies.

`publication/plot_data_generation` contains the scripts that create the
`*_plot_data.mat` files used by the plots.

In the generation scripts:

- `USE_SEEDED_NETWORKS = true` uses the seeded networks. The publication
  generators default to `SEED = 0`.
- `USE_SEEDED_NETWORKS = false` uses older single-network artifacts in
  `trained_networks`, if present.
- Dynamical-system publication data should use the active derivative-field
  seeded networks in `trained_networks/seeded_state_random_derivative_ode45_lr_0p01`,
  or a matching continuation network set if selected by the generator.
- `generate_dynamical_system_plot_data.mlx` now uses the same derivative-field
  ODE45 closed-loop method as the active DS test script: the true trajectory is
  generated from the analytic vector field, and the learned trajectory is
  generated from `dx_norm/dt = network(x_norm) - x_norm`.
- For dynamical-system plot data, `USE_CONTINUATION_IF_AVAILABLE = true` lets
  the generator use a continuation-trained network for the selected seed when
  one exists. `CONTINUATION_NETWORK_SETS = "auto"` scans continuation folders
  and prefers the largest extra-epoch continuation.
- `CLOSED_LOOP_OUTPUT_DT = []` leaves ode45 on adaptive output rows. Set it to a
  positive scalar only when a fixed plotting output grid is intentionally needed.
- `FORCE_REGENERATE = false` keeps an existing plot-data MAT file.
- `FORCE_REGENERATE = true` overwrites the plot-data MAT file for that task.
- For dynamical-system plot data, `MANUAL_INITIAL_CONDITION = []` uses the
  default IC from `dynamics_list.xlsx`; setting a vector such as `[0.01 0 0]`
  uses that base IC, with each system taking the first entries it needs.

Run individual generators when you only need one task:

- `generate_mnist_plot_data.mlx`
- `generate_toyota_plot_data.mlx`
- `generate_dynamical_system_plot_data.mlx`
- `generate_pong_plot_data.mlx`
- `generate_lqr_plot_data.mlx`
- `generate_bias_histogram_plot_data.mlx`

`generate_all_plot_data.mlx` attempts all available plot-data tasks.
