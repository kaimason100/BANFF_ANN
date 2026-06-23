# Recurrent Dale/Low-Rank RNN Reference Code

This folder contains recurrent reference code used for comparison with the
derivative-field dynamical-system workflow in the main package. The files in
this folder are preserved as reference material and are not part of the active
seeded feedforward-rate-network training pipeline.

This folder is the recurrent half of the top-level repository split. The active
feedforward package is in `../../Feedforward_network`.

The active package code does not edit these files and does not require users to
change them. They are included so the derivative-field training convention can
be inspected alongside the original RNN-style implementation.

## Files

- `rnn_learning.m` - reference training script for the low-rank/Dale-style RNN
  implementation.
- `bias_gradient_descent.m` - bias-update helper used by the reference code.
- `net_config.m` - network configuration helper.
- `net_int.m` and `net_int_large.m` - network integration wrappers.
- `int_dyn.m` - dynamical-system vector-field definitions.
- `dynamics_list.m` and `dynamics_list.xlsx` - task definitions and initial
  conditions used by the reference scripts.
- `readme.docx` - original reference-code notes.

## Relationship To The Main Package

The active dynamical-system implementation in this repository is the seeded
state-random derivative-field workflow:

```matlab
network(x_norm) = x_norm + dx_norm/dt
dx_norm/dt = network(x_norm) - x_norm
```

Training uses fixed random samples in normalized state space and a fixed
derivative-field validation set. Final testing integrates the learned vector
field with `ode45`. This was aligned with the principle of the reference RNN
code while retaining this package's feedforward architecture, deterministic
custom initializers, seeded train/validation samples, and bias-only training.

## Use

To run the recurrent reference code directly, open its scripts from this folder
and follow `readme.docx`. The feedforward release scripts do not call
`rnn_learning.m` during standard training, testing, or publication plot
generation.

## Reproducibility Note

This recurrent code is documented as-is. It has not been refactored into the
main package style because changing it would make it less useful as a reference
for comparison.
