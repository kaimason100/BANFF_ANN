# Project Layout

This repository contains two MATLAB code packages:

- `Feedforward_network` contains the seeded non-spiking feedforward/rate-network
  package, including training scripts, tests, publication plotting scripts,
  source helpers, and documentation.
- `Recurrent_network` contains the low-rank/Dale RNN reference implementation
  used for comparison with the derivative-field dynamical-system workflow.

The code files themselves were moved without changing their MATLAB
implementation. Start with `Feedforward_network/README.md` for the active
feedforward package and `Recurrent_network/README.md` for the recurrent code.

## GitHub Scope

The GitHub release includes the active code and documentation in
`Feedforward_network` and `Recurrent_network`.

Generated outputs, trained networks, local dataset MAT files, publication
plot-data MAT files, videos, and zip files are also ignored so
the repository remains code-focused and reproducible from source.

## Partner Publication

This software accompanies a partner bioRxiv preprint:

[bioRxiv 10.1101/2025.10.05.680523v1](https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1.abstract)

Please cite the linked preprint when using this code for publication-related
analyses, alongside any repository citation information provided by the authors.
