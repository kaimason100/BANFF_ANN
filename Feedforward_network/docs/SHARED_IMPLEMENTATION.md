# Shared Implementation Branch

This branch changes only the feedforward package. The recurrent package is
unchanged. It is not a new publication release and does not replace the archived
publication results. Local MATLAB and Slurm jobs now call the **same functions**;
there is no separate numerical ARC implementation to synchronise.

## Start Here

Use MATLAB R2023a or later, Deep Learning Toolbox, and Statistics and Machine
Learning Toolbox. The arm task also needs Control System Toolbox. Image
preparation may need Image Processing Toolbox. A supported GPU and Parallel
Computing Toolbox are optional. Prepared datasets are supplied in `data`;
see [DATASETS.md](DATASETS.md) for sources, fields and imputation details.

Set MATLAB's current folder to `Feedforward_network`, then run:

```matlab
run('setup.m')
run('examples/train_publication.mlx')
```

This Live Script is the single entry point for the full publication schedule.
Its `TASKS` and `SEEDS` settings restrict the run without changing a selected
task/seed's schedule. An empty `TASKS` selects all 24 tasks; `SEEDS = 0:9`
selects ten seeds. Local live plots are enabled by default. This is substantial
training, not a smoke test. Existing output files for the same task/seed are
overwritten by a fresh run; copy valuable results elsewhere first.

The schedule is defined in `src/+banff/publicationSchedule.m`:

| Tasks | Initial epochs | Learning rate |
| --- | ---: | ---: |
| Tabular classification, regression and seven MNIST-family tasks | 1,000 | 0.01 |
| Pong | 20,000 | 0.001 |
| Two-joint motor control | 10,000 | 0.01 |
| Nine derivative-field dynamical systems | 100,000 | 0.01 |

The initial DS learning rate is **0.005** for MO5 seed 3, MO13 seed 3,
MO7 seed 8, and Rikitake seeds 0 and 4. Every MO0 seed then receives a
50,000-epoch continuation at 0.001, automatically, after its initial stage.
The continuation starts from the saved best-validation network, not the last
epoch, and starts fresh Adam moments, as in the existing continuation workflow.
There are 240 initial stages and ten continuation stages.

To run a selected publication combination without opening the Live Script:

```matlab
banff.runPublication("MO0",4,true)
```

Individual task Live Scripts remain in `examples`. These are convenient for
experiments with explicit options; use `train_publication.mlx` when the exact
publication schedule, including DS exceptions, is required. Reusable functions
are `.m` files in `src`; paired `.m` copies of Live Scripts have been removed.

## Slurm

Transfer this **entire feedforward package**, including `src`, Live Scripts,
spreadsheet and prepared `data`, to ARC. Submit from `Feedforward_network`.
The shell file does not need execute permission when submitted with `sbatch`.

All 24 tasks and ten seeds, including the MO0 continuation in the same job:

```bash
sbatch --array=0-239 --export=ALL,BANFF_MODE=publication --partition=gpu-h100,gpu-a100,gpu-l40,gpu-v100 slurm/train.slurm
```

One system and all seeds, with its complete publication schedule:

```bash
sbatch --export=ALL,BANFF_MODE=publication,BANFF_TASK=MO0 slurm/train.slurm
```

One chosen seed: add `--array=4`. Clear stale `BANFF_TASK`, `BANFF_SEED`,
`BANFF_LR` and `BANFF_EPOCHS` environment variables before an all-task submission.
Publication mode rejects learning-rate/epoch overrides. The default requested
resources are one GPU, eight CPU cores, 32 GiB host RAM and a 24-hour wall-time
ceiling **per array element**. These are conservative starting requests for
the 16,000-neuron dense network, not measured minimum requirements. Host RAM is
separate from GPU memory. Inspect `sacct` elapsed time and MaxRSS after a pilot
before reducing requests; do not assume every task finishes within 24 hours.
There is no automatic full-state checkpoint/requeue mechanism. DS loss-history
checkpoints cannot resume the optimiser. Adjust resources to the cluster's
current limits (`arc.limits`). A partition list requests one allocation per element, not four
duplicate runs; Slurm policy decides scheduling and does not guarantee strict
left-to-right preference. `MATLAB_MODULE` can override `matlab/r2023a`.

For a deliberate experimental initial run, use `BANFF_MODE=train`, `BANFF_TASK`,
and optional `BANFF_LR`/`BANFF_EPOCHS`. For DS continuation use
`BANFF_MODE=continue`; its defaults are 50,000 additional epochs and 0.001.
Do not submit two jobs writing the same task/seed/output folder simultaneously.

## Saved Networks And Testing

Classification, regression, images and control now save to
`trained_networks/seeded_grouped/<task>/<task>_seed_000_network.mat`.
Control task names are `Pong` and `LQR_two_link_arm`. New tests intentionally
do not silently use legacy models without the required audit metadata.

DS retains `seeded_state_random_derivative_ode45_lr_0p01` for the base set,
including the explicitly recorded 0.005 exceptions. Continuation is saved
separately in `seeded_state_random_derivative_ode45_continued_lr_0p001_extra_50000`.
The DS test can prefer a continuation model or force the base model.

After training, run the required task Live Scripts in `tests`, or:

```matlab
run('tests/run_all_seeded_tests.mlx')
```

Control and DS tests are closed-loop: the network drives its own environment
or learned vector field. Ground truth is used for scoring, not fed back as the
network's next state. Publication data generators run after training/testing,
then their corresponding scripts in `publication/plots`. Keep
`USE_SEEDED_NETWORKS = true`; only seeded models are supported on this branch.
Regenerate old cached plot data (`FORCE_REGENERATE = true`) before comparing
new models. Default representative publication seed remains 9. Plot styling
and firing-rate colour limits are not intentionally changed.

## Leakage Checks And Deliberate Differences

* Tabular data: remove repeated complete predictor/target records, then split
  by identical predictors. Conflicting targets for the same predictors stay
  together. Nominal train/validation/test proportions remain 64/16/20, by
  predictor group. Retained/removed source indices are saved.
* Toyota has 39 duplicate complete records and 114 repeated predictor rows in
  the supplied 6,738 rows. Iris has one duplicate complete record. Their
  training/splits change deliberately; old published metrics must not be
  presented as measurements of these corrected models.
* Images: keep identical training images in one partition, remove complete
  labelled repeats, and reject overlap with the official test set. Repeated
  complete official-test records are counted only once.
* All dataset partitions are checked again after normalisation and conversion
  to single precision. Invalid indices, cross-partition identical inputs,
  altered dataset SHA-256 hashes, and mismatched training-only statistics fail
  with errors, not warnings or silent fallbacks.
* Motor episodes are partitioned by episode; training, validation and test
  targets must be distinct. Saved raw training data allows normalisation to
  be checked at test time. Effective training/validation input overlap fails.
* Pong saves its supervised sample audit and checks its train/validation
  inputs. Its final test is a separate closed-loop game, not evaluation on
  supervised training labels.
* DS uses fixed random training/validation states generated by separate
  streams. Effective overlap is an error. Continuation must retain the source
  task/seed, sample seeds, counts, range and scaling. Tests check that the new
  initial condition is not a supervised training or validation state.

Revisiting a physical state during an autonomous rollout is not itself data
leakage. Test metrics must not be used to choose parameters or repeatedly tune
a supposedly held-out test. These checks cannot detect two different records
belonging to the same real-world entity without identifiers, or prevent manual
test-driven model selection. They are not a guarantee against every possible
scientific source of leakage.

## Reproducibility And Verification

The original frozen-weight initialisers, optimiser formulae, preprocessing
formulae and seed resets are retained. No-duplicate tabular/image splits keep
the original partition calls and order. Both platforms execute the shared
functions; only plotting visibility differs in headless use. Cross-platform
bitwise training or chaotic ODE trajectory identity is **not guaranteed**:
GPU/CPU kernels, MATLAB versions, floating-point summation and adaptive solver
decisions can differ. Grouping or a failed hard check can intentionally require
different training data or prevent training.

`tests/test_shared_invariants.m` contains focused MATLAB unit tests for grouping,
duplicate failures, unchanged duplicate-free split/RNG behaviour and the full
schedule. These tests and full training have **not been executed** during this
refactor because MATLAB execution was not authorised. Static source parsing,
Live Script archive/XML checks and shell syntax checks are narrower checks,
not publication-result reproduction.

## Batch Size Metadata

For tabular tasks, an empty `MiniBatchSize` selects the historical
dataset-specific value (Iris 150; breast cancer and car quality 256;
mushroom 4096; Abalone the loaded row count; Toyota the smaller of 4096 and
the loaded row count). The resolved value is saved in
`metadata.trainingOptions.MiniBatchSize`. An explicit positive integer
override is honoured. Existing saved metadata are not rewritten.
