# BANFF Neural-Network Experiments

MATLAB code accompanying *Rapidly Reconfigurable Dynamic Computing in Neural
Networks with Fixed Synaptic Connectivity*. BANFF networks retain fixed
synaptic weights and learn task-specific neuronal biases.

The repository has two independent packages:

- `Feedforward_network` contains the seeded feedforward experiments for
  classification, regression, autonomous dynamical systems, Pong, and
  two-joint motor control.
- `Recurrent_network` contains the low-rank recurrent experiments, including
  multisystem learning, bias switching, representational drift, and the
  Dale-constrained network.

Start with the README inside the package you intend to run. Feedforward users
should then follow `Feedforward_network/docs/GETTING_STARTED.md` and
`Feedforward_network/docs/DATASETS.md` in order.

## Publication

This repository contains the software accompanying the article:

> Kai Mason, Sonia Sennik, Claudia Clopath, Aaron Gruber, and Wilten Nicola.
> *Rapidly Reconfigurable Dynamic Computing in Neural Networks with Fixed
> Synaptic Connectivity*. **Communications Biology** (accepted; final volume,
> article number, and journal DOI pending).

Wilten Nicola is the corresponding author. Until the final journal record is
available, cite the [bioRxiv preprint](https://doi.org/10.1101/2025.10.05.680523)
alongside this repository. `CITATION.cff` provides machine-readable citation
metadata for GitHub and reference managers.

## Repository Scope

GitHub contains source code and documentation. Third-party datasets must be
obtained from their upstream sources and prepared locally. Trained networks,
generated results, publication plot-data MAT files, figures, videos, and local
compute-environment packages are intentionally excluded. Dataset sources,
preparation details, and third-party licence terms are documented in
`Feedforward_network/docs/DATASETS.md` and `THIRD_PARTY_NOTICES.md`.

## Reproducibility

All explicit random seeds, dataset splits, and train-fitted preprocessing
statistics are fixed. Exact bitwise equality across operating systems, MATLAB
versions, CPUs, and GPUs is not guaranteed because floating-point kernels and
adaptive ODE integration may differ. Long chaotic rollouts can amplify those
small differences. Compare saved parameters directly when checking
initialisation, and use aggregate task metrics from a controlled environment
for publication comparisons.

## Requirements

The release targets MATLAB R2023a or later. Feedforward experiments require the
Deep Learning Toolbox and Statistics and Machine Learning Toolbox. MNIST-family
image preparation uses Image Processing Toolbox functions, and the two-joint
arm task requires Control System Toolbox. A supported GPU is optional. Recurrent
experiments use standard MATLAB numerical solvers and the supplied source files.

## License

Except for identified third-party material, this software is released under
the MIT License. See `LICENSE`. Third-party attribution and applicable licence
terms, including those for the datasets used by the experiments and `slanCM`
plotting compatibility support, are provided in `THIRD_PARTY_NOTICES.md` and
`LICENSES`.
