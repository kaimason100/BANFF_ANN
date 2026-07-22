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

This repository contains the software accompanying the journal publication
*Rapidly Reconfigurable Dynamic Computing in Neural Networks with Fixed
Synaptic Connectivity*. Please cite the journal article using its final
bibliographic record.

## Repository Scope

GitHub contains source code and documentation. Trained networks, source dataset
MAT files, generated results, publication plot-data MAT files, figures, videos,
and local compute-environment packages are intentionally excluded. These files
must be generated or supplied locally as described in the package documentation.

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
