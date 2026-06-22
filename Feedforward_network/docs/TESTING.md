# Testing

For the full end-to-end run order, start with `docs/GETTING_STARTED.md`. This
file focuses only on saved-network tests.

Run active seeded tests from the repository root, `tests`, or any subfolder:

```matlab
run("tests/run_all_seeded_tests.mlx")
```

The active tests load saved seeded artifacts from `trained_networks`; they do
not retrain networks.

## Active aggregate tests

- `test_classification_seeded.mlx`
- `test_regression_seeded.mlx`
- `test_mnist_family_seeded.mlx`
- `test_dynamical_systems_seeded_state_random_derivative_ode45.mlx`
- `test_motor_control_seeded.mlx`
- `test_pong_seeded.mlx`

Individual per-task tests are also available for tabular, regression, and
MNIST-family tasks.

## Hard checks

The tests assert the relevant saved metadata before reporting metrics:

- seed value and seed list `0:9`
- disjoint train/validation/test splits where applicable
- train-fitted normalization statistics
- official MNIST-family test split usage
- motor-control train/test target separation
- Pong closed-loop simulation settings
- dynamical-system derivative-field metadata matching the active local trainer

Regression tests quote RMSE, Pearson `r`, and Pearson p-value as mean and SD
across the saved seeds. Classification and MNIST-family tests report accuracy
across seeds. Closed-loop task tests report task-appropriate closed-loop
metrics and seed summaries.

## Dynamical systems

The active DS test is `test_dynamical_systems_seeded_state_random_derivative_ode45.mlx`.
It expects networks saved by the derivative-field trainer and tests them by
integrating the learned continuous-time vector field with `ode45`:

- network equation `dx_norm/dt = network(x_norm) - x_norm`
- training supervisor `network(x_norm) = x_norm + dx_norm/dt`
- total epochs `100000`
- validation interval `200`
- fixed training samples `1000`
- fixed validation samples `100`
- sample RNG seed `12345`
- validation RNG seed `9000`
- reference scaling sample time `1000`
- default closed-loop rollout length `1000` time units
- default `ode45` output grid `0.01` time units; this is an output grid, not a
  fixed-step integrator

Publication DS networks were initially trained for `100000` epochs. Most
task/seed combinations used learning rate `0.01`; the exceptions were `MO5`
seed `3`, `MO13` seed `3`, `MO7` seed `8`, and `Rikitake` seeds `0` and `4`,
which used learning rate `0.005`. All `MO0` seeds were then continued for
`50000` extra epochs at learning rate `0.001`.

By default the test automatically uses a matching continuation network if one
exists for that task/seed, then falls back to the base derivative-field network.
The results table records `NetworkSet`, `UsedContinuation`, and `NetworkPath`.

## Numerical reproducibility

Tests are CPU-only where practical and all explicit random seeds are fixed.
Small numerical differences can still occur across machines because MATLAB may
use different CPU math libraries and ODE step decisions on different platforms.
For chaotic closed-loop dynamical systems, those differences can grow over long
rollouts. This is expected numerical sensitivity rather than a change in saved
weights or data splits.

## Weight consistency

`check_seeded_network_weight_consistency.mlx` verifies same-seed frozen-weight
consistency across saved seeded networks. Hidden-to-hidden weights must match
exactly; encoder and decoder weights are compared over the largest shared
leading submatrix. The active checker includes the main seeded networks, the
derivative-field dynamical-system networks, and any derivative continuation
folders present locally. It uses a dominant largest-matrix reference for
encoder/decoder comparisons when possible, with all-pairs fallback only for
unusual non-nested matrix shapes.
