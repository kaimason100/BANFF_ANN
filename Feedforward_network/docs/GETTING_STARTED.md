# Getting Started With The Feedforward Package

This guide is written for a new user who has not seen the code before. It
explains what the package contains, how to prepare MATLAB, how to run the
training scripts, how to test saved networks, and how to regenerate publication
plot data.

Run the feedforward scripts from inside `Feedforward_network`.

This package accompanies the partner bioRxiv preprint
https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1.abstract. The
GitHub repository is intended to contain source code and documentation; large
local datasets, trained networks, and generated outputs are
excluded from normal commits.

## 1. Software requirements

Use MATLAB with the Deep Learning Toolbox and Statistics and Machine Learning
Toolbox. Newer MATLAB versions should run the local scripts, but small numerical
differences can occur across MATLAB versions, operating systems, CPUs, and GPUs.

The scripts are MATLAB live scripts (`.mlx`). They can be opened in Live Editor
or run from the command window with `run(...)`.

Recommended local workflow:

```matlab
cd "path/to/GitHub"
run("tests/run_all_seeded_tests.mlx")
```

The scripts try to locate the repository root automatically, so they can usually
be run from any folder. If a publication plot-data generator says it cannot find
the repository root, set `MANUAL_REPO_ROOT` at the top of that generator to the
absolute path of this repository.

## 2. Folder map

- `examples/classification` contains seeded classification training scripts.
  This includes tabular tasks and all MNIST-family image tasks.
- `examples/regression` contains seeded regression training scripts.
- `examples/control` contains the seeded two-joint LQR motor-control trainer.
- `examples/pong` contains the seeded Pong trainer and the playback script.
- `examples/dynamical_systems` contains the active seeded derivative-field
  dynamical-system trainers and `dynamics_list.xlsx`.
- `tests` contains saved-network test scripts. These load saved networks and
  never retrain.
- `src/initializers` contains the deterministic frozen-weight initializers.
- `src/preprocessing` contains reusable time-delay helpers.
- `data` contains the datasets used by local training and test scripts.
- `trained_networks` is the expected output folder for saved networks.
- `publication/plot_data_generation` creates the MAT files used by publication
  plots.
- `publication/plots` contains the legacy publication figure scripts.

## 3. First-time setup

1. Put the repository somewhere MATLAB can access.
2. Put required dataset files in `data`. See `docs/DATASETS.md` for exact file
   names, variable names, and expected array shapes.
3. Open MATLAB and set the current folder to the repository root.
4. Run a quick path/root check:

```matlab
run("tests/check_seeded_network_weight_consistency.mlx")
```

This check requires saved networks. If you have not trained anything yet, it may
fail because the artifacts are missing. That is expected before training.

## 4. Training order

All active training scripts use the seeded convention `SEEDS = 0:9`. For a given
seed, the frozen physical weight matrices are initialized consistently across
tasks up to the largest shared input/output submatrix. Training learns biases
while the custom weight matrices are frozen.

The safest full local order is:

1. Prepare datasets in `data`.
2. Train tabular classification tasks.
3. Train MNIST-family image tasks.
4. Train regression tasks.
5. Train Pong.
6. Train motor control.
7. Train dynamical systems.
8. Run tests.
9. Generate publication plot data.
10. Run publication plotting scripts.

You can run individual tasks instead of the batch runners when you only need one
network family.

## 5. Classification training

Tabular classification tasks:

```matlab
run("examples/classification/train_classification_seeded.mlx")
```

This batch runner trains:

- `iris`
- `breast_cancer`
- `car_quality`
- `mushroom`

Individual scripts are also available:

```matlab
run("examples/classification/train_iris_seeded.mlx")
run("examples/classification/train_breast_cancer_seeded.mlx")
run("examples/classification/train_car_quality_seeded.mlx")
run("examples/classification/train_mushroom_seeded.mlx")
```

MNIST-family image tasks:

```matlab
run("examples/classification/train_mnist_family_seeded.mlx")
```

This batch runner trains:

- `mnist`
- `mnist_fashion`
- `kmnist`
- `afro_mnist_ethiopic`
- `afro_mnist_nko`
- `afro_mnist_osmanya`
- `afro_mnist_vai`

Individual scripts are also available:

```matlab
run("examples/classification/train_mnist_seeded.mlx")
run("examples/classification/train_mnist_fashion_seeded.mlx")
run("examples/classification/train_kmnist_seeded.mlx")
run("examples/classification/train_afro_mnist_ethiopic_seeded.mlx")
run("examples/classification/train_afro_mnist_nko_seeded.mlx")
run("examples/classification/train_afro_mnist_osmanya_seeded.mlx")
run("examples/classification/train_afro_mnist_vai_seeded.mlx")
```

Outputs are saved to:

```text
trained_networks/seeded/<task>/<task>_seed_<seed>_network.mat
```

## 6. Regression training

Run both regression tasks:

```matlab
run("examples/regression/train_regression_seeded.mlx")
```

Or run one task:

```matlab
run("examples/regression/train_abalone_seeded.mlx")
run("examples/regression/train_toyota_seeded.mlx")
```

Regression tests report only aggregate RMSE, Pearson `r`, and Pearson p-value
as mean and SD across seeds.

Outputs are saved to:

```text
trained_networks/seeded/<task>/<task>_seed_<seed>_network.mat
```

## 7. Pong training

Run:

```matlab
run("examples/pong/train_pong_seeded.mlx")
```

The Pong trainer generates its own supervised training samples inside the
script, then saves seeded networks to:

```text
trained_networks/seeded/pong/pong_seed_<seed>_network.mat
```

The Pong test is closed-loop: it simulates games and counts successful network
paddle contacts and misses. To play or visualize Pong separately, use:

```matlab
run("examples/pong/playPong.mlx")
```

## 8. Motor-control training

Run:

```matlab
run("examples/control/train_motor_control_seeded.mlx")
```

The motor-control script generates supervised LQR episodes and saves seeded
networks to:

```text
trained_networks/seeded/motor_control/motor_control_seed_<seed>_network.mat
```

The test is closed-loop: each saved network controls the arm toward held-out
targets rather than being evaluated only on precomputed open-loop samples.

## 9. Dynamical-system training

The active dynamical-system workflow is the seeded derivative-field ODE45
implementation. It trains on fixed random normalized state-space samples and
uses the true dynamical-system equations to build the supervisor:

```matlab
network(x_norm) = x_norm + dx_norm/dt
```

Testing integrates the learned field:

```matlab
dx_norm/dt = network(x_norm) - x_norm
```

Run the full local batch:

```matlab
run("examples/dynamical_systems/train_dynamical_systems_seeded_state_random_derivative_ode45.mlx")
```

Run one specific task/seed:

```matlab
run("examples/dynamical_systems/train_dynamical_system_single_seed_state_random_derivative_ode45.mlx")
```

At the top of that script, edit:

```matlab
TASK_NAME = "Lorenz";
SEED = 0;
TOTAL_EPOCHS = 100000;
LEARNING_RATE = 0.01;
```

Use the single-task script for targeted reruns where the publication schedule
does not use the default learning rate. If you are recreating the publication
artifact tree exactly, also set `OUTPUT_NETWORK_SET` to the network set you want
the test scripts to load. The current publication DS artifacts use the active
base set name `seeded_state_random_derivative_ode45_lr_0p01`; the saved metadata
inside each `.mat` file records the actual learning rate used for that
task/seed.

Publication initial-training schedule:

| Task | Seeds | Initial epochs | Initial learning rate |
| --- | --- | ---: | ---: |
| Lorenz | `0:9` | `100000` | `0.01` |
| MO0 | `0:9` | `100000` | `0.01` |
| MO5 | all except `3` | `100000` | `0.01` |
| MO5 | `3` | `100000` | `0.005` |
| MO7 | all except `8` | `100000` | `0.01` |
| MO7 | `8` | `100000` | `0.005` |
| MO13 | all except `3` | `100000` | `0.01` |
| MO13 | `3` | `100000` | `0.005` |
| Rikitake | all except `0, 4` | `100000` | `0.01` |
| Rikitake | `0, 4` | `100000` | `0.005` |
| SprottB | `0:9` | `100000` | `0.01` |
| SprottC | `0:9` | `100000` | `0.01` |
| SprottS | `0:9` | `100000` | `0.01` |

Continue selected trained networks with a lower learning rate:

```matlab
run("examples/dynamical_systems/continue_dynamical_systems_seeded_state_random_derivative_ode45.mlx")
```

At the top of the continuation script, edit `continuationPairs`, for example:

```matlab
continuationPairs = {
    'MO0', 0
    'MO0', 1
    'MO0', 2
    'MO0', 3
    'MO0', 4
    'MO0', 5
    'MO0', 6
    'MO0', 7
    'MO0', 8
    'MO0', 9
    };
```

The publication continuation stage was run for all `MO0` seeds only: `50000`
extra epochs with learning rate `0.001`. The local continuation default is
`50000` extra epochs unless you edit `CONTINUATION_EPOCHS`.

Active task rows from `dynamics_list.xlsx`:

- `Lorenz`
- `MO0`
- `MO5`
- `MO7`
- `MO13`
- `Rikitake`
- `SprottB`
- `SprottC`
- `SprottS`

Important default parameters:

- seeds `0:9`
- total epochs `100000`
- validation every `200` epochs
- fixed training samples `1000`
- fixed validation samples `100`
- training sample RNG seed `12345`
- validation RNG seed `9000`
- normalized state range `[-1 1]`
- reference scaling sample time `1000`
- Adam learning rate `0.01` except for the explicit task/seed exceptions in the
  publication schedule table above
- Adam beta1 `0.99`, beta2 `0.999`, epsilon `5e-6`

Outputs are saved to:

```text
trained_networks/seeded_state_random_derivative_ode45_lr_0p01/<task>/<task>_seed_<seed>_network.mat
```

Continuation outputs are saved to:

```text
trained_networks/seeded_state_random_derivative_ode45_continued_lr_0p001_extra_<epochs>/<task>/<task>_seed_<seed>_network.mat
```

Local runs keep live plots visible where applicable. Plotting visibility does
not change the seeded training parameters.

## 10. Running tests

After training, run all active tests:

```matlab
run("tests/run_all_seeded_tests.mlx")
```

Or run aggregate tests by task family:

```matlab
run("tests/test_classification_seeded.mlx")
run("tests/test_mnist_family_seeded.mlx")
run("tests/test_regression_seeded.mlx")
run("tests/test_pong_seeded.mlx")
run("tests/test_motor_control_seeded.mlx")
run("tests/test_dynamical_systems_seeded_state_random_derivative_ode45.mlx")
```

Individual task tests are also available, for example:

```matlab
run("tests/test_mnist_seeded.mlx")
run("tests/test_abalone_seeded.mlx")
run("tests/test_toyota_seeded.mlx")
```

Tests load saved networks only. They do not retrain networks. They also enforce
hard checks for:

- seed metadata
- disjoint train/validation/test splits
- training-only normalization statistics
- correct MNIST-family official test split usage
- closed-loop settings for Pong, motor control, and dynamical systems

## 11. Dynamical-system test controls

At the top of the DS derivative-field test, you can change:

```matlab
CLOSED_LOOP_ROLLOUT_LENGTH = 1000;
ODE_OUTPUT_DT = 0.01;
TEST_IC_RANDOM_SEED = 9100;
TEST_IC_PERTURBATION_SCALE = 0;
USE_CONTINUATION_IF_AVAILABLE = true;
CONTINUATION_NETWORK_SETS = "auto";
```

`CLOSED_LOOP_ROLLOUT_LENGTH` is the physical time interval integrated by
`ode45`. `ODE_OUTPUT_DT` is the requested output grid, not a fixed-step
integrator. With `USE_CONTINUATION_IF_AVAILABLE = true`, the test scans for
matching continuation networks and uses them automatically before falling back
to the base derivative-field network.

## 12. Weight consistency checks

To check that same-seed networks share the intended frozen weights:

```matlab
run("tests/check_seeded_network_weight_consistency.mlx")
```

This verifies that hidden-to-hidden weights match exactly across tasks for a
given seed, and that encoder/decoder weights match over the largest shared
leading submatrix when tasks have different input/output sizes. The checker
loads the active seeded folder, the active derivative-field dynamical-system
folder, and any local derivative continuation folders by default. It compares
each seed group part-by-part against the largest matrix that contains the other
task matrices as leading submatrices, which avoids an unnecessary full all-pairs
comparison.

## 13. Publication plot-data generation

Publication plot scripts load MAT files from `publication/plots`. Generate or
refresh those MAT files before running the figure scripts:

```matlab
run("publication/plot_data_generation/generate_all_plot_data.mlx")
```

Or run individual generators:

```matlab
run("publication/plot_data_generation/generate_mnist_plot_data.mlx")
run("publication/plot_data_generation/generate_toyota_plot_data.mlx")
run("publication/plot_data_generation/generate_dynamical_system_plot_data.mlx")
run("publication/plot_data_generation/generate_pong_plot_data.mlx")
run("publication/plot_data_generation/generate_lqr_plot_data.mlx")
run("publication/plot_data_generation/generate_bias_histogram_plot_data.mlx")
```

Important generator settings:

- `USE_SEEDED_NETWORKS = true` loads a seeded network.
- `SEED = 0` chooses which seeded network to use.
- `USE_SEEDED_NETWORKS = false` loads older single-network artifacts, if they
  are present.
- `FORCE_REGENERATE = false` keeps existing plot-data files.
- `FORCE_REGENERATE = true` overwrites existing plot-data files.
- For dynamical systems, the plot-data generator uses the active
  derivative-field ODE45 testing method. `DYNAMICAL_SYSTEM_NETWORK_SET` selects
  the base seeded derivative network set, and
  `USE_CONTINUATION_IF_AVAILABLE = true` lets the generator use a matching
  continuation network when one exists.
- `CLOSED_LOOP_OUTPUT_DT = []` leaves dynamical-system plot-data generation on
  adaptive ode45 output rows. Set it to a positive scalar only when you
  intentionally want a fixed plotting output grid.

After generating plot data, run scripts in `publication/plots`, for example:

```matlab
run("publication/plots/plot_MNIST_confusion_matrix.mlx")
run("publication/plots/plot_toyota_correlation.mlx")
run("publication/plots/plot_lorenz_phase_portrait.mlx")
```

The plot scripts are intended to preserve legacy figure formatting.

## 14. Common errors

`Missing seeded model` means the relevant training script has not been run, or
the saved network is not in the expected `trained_networks` folder.

`Missing dataset` means a required `.mat` dataset is absent from `data`, or the
file name does not match the expected name. See `docs/DATASETS.md`.

`Could not locate repository root` usually means MATLAB Live Editor executed the
script from a temporary folder. For publication plot-data generators, set
`MANUAL_REPO_ROOT`. For other scripts, run from the repository root or from a
folder inside the repository.

Metadata assertion failures usually mean an older saved network is being tested
with a newer test script. Regenerate the network with the active seeded training
script that matches the test.

Different DS closed-loop trajectories across machines are expected for long
chaotic rollouts. The saved weights and explicit seeds can match while adaptive
ODE decisions and floating-point details still differ slightly.
