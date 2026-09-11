# Feedforward BANFF experiments

MATLAB code for the seeded feedforward experiments accompanying *Rapidly
Reconfigurable Dynamic Computing in Neural Networks with Fixed Synaptic
Connectivity*. Network weights are fixed; neuronal biases are trained.

Start with [Getting Started](docs/GETTING_STARTED.md). Dataset files, trained
networks, and generated outputs are not distributed.

## Workflow

From `Feedforward_network` in MATLAB:

```matlab
run('setup.m')
run('examples/train_publication.mlx')
run('tests/run_all_seeded_tests.mlx')
```

Training covers 24 tasks and seeds 0–9, including the dynamical-system
learning-rate exceptions and MO0 continuation. This is a full training run;
edit `TASKS` and `SEEDS` in the training entry script for a smaller selection.
Training overwrites existing models for the selected task and seed.

## Reference

- [Getting started](docs/GETTING_STARTED.md): setup, training, testing and figures.
- [Datasets](docs/DATASETS.md): sources, formats, checksums and licences.
- [Testing](docs/TESTING.md): evaluation settings and integrity checks.
- [Publication figures](publication/README.md): plot-data generation and plotting.
- [Package structure](docs/PACKAGE_STRUCTURE.md): source and entry-point locations.

Exact numerical results can vary with MATLAB and hardware, particularly for
chaotic closed-loop rollouts. Retain the evaluated models and test outputs
when recording publication results.

The repository MIT licence covers the code except identified third-party
material. Dataset and plotting-library terms are listed in
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
