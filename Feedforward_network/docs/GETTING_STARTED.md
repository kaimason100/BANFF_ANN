# Getting started

Run these commands from the `Feedforward_network` folder in MATLAB.
Dataset files, models, figures, and generated results are produced or prepared
locally.

## Setup

Use MATLAB R2023a or later with Deep Learning Toolbox and Statistics and
Machine Learning Toolbox. The arm task requires Control System Toolbox.
MNIST-family image preparation uses Image Processing Toolbox. GPU execution
is optional.

```matlab
run('setup.m')
```

Download sources and prepare the required MAT files by following
[DATASETS.md](DATASETS.md). Do not pre-normalize the prepared data: training
fits normalization statistics on training samples and saves them for testing.

## Train

```matlab
run('examples/train_publication.mlx')
```

This runs the complete 24-task, ten-seed publication schedule, including
DS learning-rate exceptions and the MO0 continuation. Edit `TASKS` and `SEEDS`
in that Live Script to select a subset. Existing files for those task/seed
pairs are overwritten. This is substantial training, not a quick test.

For a selected task and seed with the same publication schedule:

```matlab
banff.runPublication("Pong",9,true)
```

Task-specific Live Scripts are also available under `examples/classification`,
`examples/regression`, `examples/control`, `examples/pong` and
`examples/dynamical_systems`. Use the publication entry point when the exact
schedule is required.

Dataset and control models are saved under `trained_networks/seeded_grouped`.
DS base models use `seeded_state_random_derivative_ode45_lr_0p01`; MO0
continuations use `seeded_state_random_derivative_ode45_continued_lr_0p001_extra_50000`.

## Test saved models

After training:

```matlab
run('tests/run_all_seeded_tests.mlx')
```

This evaluates all default tasks and seeds. For a partial training run, use
the corresponding individual tests instead:

```matlab
banff.test('iris',9)
run('tests/test_motor_control_seeded.mlx')
run('tests/test_pong_seeded.mlx')
run('tests/test_dynamical_systems_seeded_state_random_derivative_ode45.mlx')
```

Edit `SEEDS` in the motor, Pong or DS test to select models; the DS test also
has `TASKS`. Tests load saved models and do not retrain. Motor, Pong and DS
tests are closed-loop. DS testing prefers an available matching continuation
unless `USE_CONTINUATION_IF_AVAILABLE` is disabled.

[TESTING.md](TESTING.md) documents metrics, ICs, continuation selection and
numerical settings. Missing models or failed integrity checks stop evaluation.
The optional `tests/check_seeded_network_weight_consistency.mlx` compares fixed
weights across saved tasks. Fast implementation checks can be run separately:

```matlab
results = runtests('tests/test_shared_invariants.m');
assertSuccess(results)
```

## Generate publication figures

```matlab
run('publication/plot_data_generation/generate_all_plot_data.mlx')
```

Generators default to seed 9. Change `SEED` deliberately and
regenerate cached files after changing models or settings. Individual
generators are available for MNIST, Toyota, DS, Pong, the arm and bias histograms.
Then run the required figure script under `publication/plots`.
See [publication/README.md](../publication/README.md) for the controls.

If automatic root discovery fails, set `MANUAL_REPO_ROOT` in the generator to
the absolute package path. A fresh checkout has no trained networks or plot
data, so those stages must be completed before plotting.

## Reproducibility

Keep the source revision, evaluated model files, test settings and outputs
together. Fixed random seeds do not guarantee bitwise equality across MATLAB
versions or hardware. Chaotic DS trajectories can amplify small numerical
differences; compare aggregate metrics under a controlled environment.
