# Feedforward Non-spiking Neural Networks Bias Learning

MATLAB live-script release for seeded non-spiking/rate-network experiments. The
active package now contains only the seeded workflows used for the publication
runs.

## Partner Publication

This feedforward package accompanies the partner bioRxiv preprint:

[bioRxiv 10.1101/2025.10.05.680523v1](https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1.abstract)

Use the scripts and documentation here to reproduce the seeded feedforward
experiments associated with that manuscript. Exact trained-network artifacts and
large generated outputs are intentionally not part of the GitHub code release.

## Start here

New users should read these two documents first:

- `docs/GETTING_STARTED.md` - setup, run order, training, testing,
  publication plotting, and common errors.
- `docs/DATASETS.md` - dataset sources and the exact `.mat` formats
  expected by the scripts.

The shorter sections below summarize the active package.

## Repository layout

- `examples/classification` - seeded tabular and MNIST-family classification
  training scripts.
- `examples/regression` - seeded regression training scripts.
- `examples/dynamical_systems` - seeded derivative-field dynamical-system
  training, single-network reruns, and continuation training.
- `examples/control` - seeded LQR two-link-arm motor-control training.
- `examples/pong` - seeded Pong training plus playback.
- `tests` - seeded saved-network tests; these load saved artifacts and do not
  retrain.
- `src/initializers` - deterministic custom weight initializers.
- `src/preprocessing` - reusable time-delay helpers.
- `src/dynamical_systems` - dynamical-system right-hand-side helpers.
- `trained_networks` - saved seeded network artifacts.
- `data` - datasets expected by local scripts.
- `publication` - publication plot scripts and plot-data generation scripts.
- `docs` - newcomer instructions, dataset preparation, package notes, and
  testing notes.

## Active training scripts

All active training scripts use the seeded release convention `SEEDS = 0:9` and
save to `trained_networks/seeded/<task>` unless stated otherwise.

Classification:

- `train_iris_seeded.mlx`
- `train_breast_cancer_seeded.mlx`
- `train_car_quality_seeded.mlx`
- `train_mushroom_seeded.mlx`
- `train_classification_seeded.mlx` batch runner

These tabular classification scripts train for `1000` epochs with Adam
learning rate `0.01`.

MNIST-family classification, also under `examples/classification`:

- `train_mnist_seeded.mlx`
- `train_mnist_fashion_seeded.mlx`
- `train_kmnist_seeded.mlx`
- `train_afro_mnist_ethiopic_seeded.mlx`
- `train_afro_mnist_nko_seeded.mlx`
- `train_afro_mnist_osmanya_seeded.mlx`
- `train_afro_mnist_vai_seeded.mlx`
- `train_mnist_family_seeded.mlx` batch runner

These MNIST-family scripts train for `1000` epochs with Adam learning rate
`0.01`.

Regression:

- `train_abalone_seeded.mlx`
- `train_toyota_seeded.mlx`
- `train_regression_seeded.mlx` batch runner

These regression scripts train for `1000` epochs with Adam learning rate
`0.01`.

Closed-loop tasks:

- `train_motor_control_seeded.mlx`
- `train_pong_seeded.mlx`
- `train_dynamical_systems_seeded_state_random_derivative_ode45.mlx`
- `train_dynamical_system_single_seed_state_random_derivative_ode45.mlx`
- `continue_dynamical_systems_seeded_state_random_derivative_ode45.mlx`

Pong trains for `20000` epochs with Adam learning rate `0.001`. Motor control
trains for `10000` epochs with Adam learning rate `0.01`.

The active dynamical-system derivative-field trainer saves to
`trained_networks/seeded_state_random_derivative_ode45_lr_0p01/<task>` by
default. It trains the network on fixed random state-space samples with
supervisor `network(x_norm) = x_norm + dx_norm/dt`, selects parameters by fixed
validation derivative-field MSE, and tests by integrating
`dx_norm/dt = network(x_norm) - x_norm` with `ode45`.

- task rows `[5, 34, 39, 41, 47, 32, 8, 9, 25]`
- seeds `0:9`
- total epochs `100000`
- validation every `200` epochs
- fixed validation samples `100`
- validation RNG seed `9000`
- reference scaling sample time `1000`
- fixed training samples `1000`
- training sample RNG seed `12345`
- normalized state range `[-1 1]`
- Adam beta1 `0.99`, beta2 `0.999`, epsilon `5e-6`

Publication DS learning-rate schedule:

- Initial derivative-field training was `100000` epochs for every task/seed.
- The initial learning rate was `0.01` for all task/seed combinations except:
  `MO5` seed `3`, `MO13` seed `3`, `MO7` seed `8`, and `Rikitake` seeds `0`
  and `4`, which used learning rate `0.005`.
- All `MO0` seeds `0:9` were then continued for an extra `50000` epochs with
  learning rate `0.001`. The local continuation script keeps the broader
  ARC-used editable continuation list so underperforming task/seed pairs can be
  reproduced without changing the training implementation.

Use `train_dynamical_systems_seeded_state_random_derivative_ode45.mlx` for the
standard `0.01` batch, `train_dynamical_system_single_seed_state_random_derivative_ode45.mlx`
for the `0.005` task/seed reruns, and
`continue_dynamical_systems_seeded_state_random_derivative_ode45.mlx` for the
MO0 continuation stage.

Continuation training loads selected derivative-field networks and continues
bias-only training with learning rate `0.001`, saving to a continuation network
set such as
`trained_networks/seeded_state_random_derivative_ode45_continued_lr_0p001_extra_50000/<task>`.
Local scripts keep live plots visible where applicable.

## Reproducibility notes

The scripts fix all explicit random sources used by the package: split seeds,
weight-seed offsets, derivative-field state-sample seeds, validation/test seeds,
and task array ordering. Local training scripts also state execution choices
explicitly where applicable, avoiding hidden MATLAB defaults. They also save
preprocessing statistics fitted on training data only and the tests assert
split/normalization metadata before reporting metrics.

Exact bitwise equality across machines is not guaranteed. MATLAB, CPU BLAS
libraries, GPU kernels, `trainNetwork`, `dlnetwork`, and adaptive `ode45` can
produce small floating-point differences across operating systems, MATLAB
releases, and hardware. Closed-loop chaotic dynamical systems amplify tiny
numerical differences over time, so trajectories and downstream metrics may
differ slightly even when the saved weights, seeds, and hyperparameters match.
For publication comparisons, use the same saved network artifacts and quote the
seeded aggregate statistics from one controlled execution environment.

## Release zip contents

Release zips are code/documentation packages only. They exclude saved networks,
generated output folders, publication plot-data MAT files, videos, local dataset
MAT files, and old release folders. New users should recreate local datasets
following `docs/DATASETS.md` and then run the active seeded scripts.

## Running tests

Run active seeded tests from the repository root or any subfolder:

```matlab
run("tests/run_all_seeded_tests.mlx")
```

Or run an individual test such as:

```matlab
run("tests/test_mnist_seeded.mlx")
run("tests/test_regression_seeded.mlx")
run("tests/test_dynamical_systems_seeded_state_random_derivative_ode45.mlx")
```

Tests load saved networks only. They do not retrain and they do not recompute
normalization from validation or test data.

## Dataset files

Place user-prepared datasets in `data` where required. See
`docs/DATASETS.md` for download locations, expected file names, expected MAT
variables, column conventions, and example conversion code. Do not pre-normalize
features before saving datasets; the training scripts fit normalization from
training data only.
