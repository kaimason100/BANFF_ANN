# BANFF_ANN

MATLAB examples for simulating **BANFF** feedforward artificial neural networks (ANNs), organised as **MATLAB Live Scripts (`.mlx`)**.

> GitHub usually won’t render `.mlx` contents (“Binary file not shown”). Download/clone the repo and run the scripts in MATLAB.

## Related manuscript

This repository is linked to the following bioRxiv preprint. **If you use this code (or derivatives of it) in academic work, please cite the manuscript:**

Mason, K., Sennik, S., Clopath, C., Gruber, A., & Nicola, W. *Rapidly Reconfigurable Dynamic Computing in Neural Networks with Fixed Synaptic Connectivity.* bioRxiv (2025). https://doi.org/10.1101/2025.10.05.680523

- Abstract: https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1.abstract  
- Full preprint page: https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1

## Citation

### APA
> Mason, K., Sennik, S., Clopath, C., Gruber, A., & Nicola, W. (2025). *Rapidly Reconfigurable Dynamic Computing in Neural Networks with Fixed Synaptic Connectivity*. bioRxiv. https://doi.org/10.1101/2025.10.05.680523

### BibTeX
```bibtex
@article{Mason2025BANFF,
  title   = {Rapidly Reconfigurable Dynamic Computing in Neural Networks with Fixed Synaptic Connectivity},
  author  = {Mason, Kai and Sennik, Sonia and Clopath, Claudia and Gruber, Aaron and Nicola, Wilten},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.10.05.680523},
  url     = {https://www.biorxiv.org/content/10.1101/2025.10.05.680523v1}
}
```

## What’s in this repo

The codebase is MATLAB-only. Most content lives under `Feedforward_Network/` and is grouped by application area (classification, regression, dynamical systems, control, and a simple game demo).

### Example notebooks (Live Scripts)

**Classification**
- `Feedforward_Network/Classification/Iris_classification.mlx`
- `Feedforward_Network/Classification/Breast_cancer_classification.mlx`
- `Feedforward_Network/Classification/Car_quality_classification.mlx`
- `Feedforward_Network/Classification/Mushroom_classification.mlx`
- `Feedforward_Network/Classification/MNIST_variants.mlx`

**Regression**
- `Feedforward_Network/Regression/Abalone_age_prediction.mlx`
- `Feedforward_Network/Regression/Car_price_regression.mlx`

**Dynamical systems / time-delay inputs**
- `Feedforward_Network/Dynamical_systems/Dynamical_systems.mlx`
- `Feedforward_Network/Dynamical_systems/extractTimeDelayedInputs.mlx`
- `Feedforward_Network/Dynamical_systems/extractTimeDelayedInputsFeedback.mlx`

**Motor control**
- `Feedforward_Network/Motor Control/motor_control.mlx`

**Game demo**
- `Feedforward_Network/Game/Pong.mlx`
- `Feedforward_Network/Game/playPong.mlx`

**Custom weight initialisation helpers**
- `Feedforward_Network/customWeights1_git.mlx`
- `Feedforward_Network/customWeights2_git.mlx`
- `Feedforward_Network/customWeights3_git.mlx`

## Requirements

- **MATLAB** (recent version recommended)
- **Deep Learning Toolbox** (for standard feedforward network functions)

Some scripts may rely on datasets being available locally or being loaded/downloaded by the script itself (check the top of each `.mlx`).

## Quick start

1. Clone or download the repository.
2. Open MATLAB.
3. Set the MATLAB **Current Folder** to the repository root (or add it to your path).
4. Open any `.mlx` file and click **Run**.

Suggested first run:
- `Feedforward_Network/Classification/Iris_classification.mlx`

## Project status

This repository is still evolving and may change structure without notice.

## Contributing

Issues and pull requests are welcome:
- If you find a bug, include MATLAB version + toolbox versions + the script you ran.
- If you add an example, place it in the relevant subfolder and keep it runnable as a standalone `.mlx`.

## Licence

No licence is currently specified in the repository. If you intend others to reuse this code, consider adding a standard open-source licence (e.g., MIT, BSD-3, GPL-3.0).

## Contact

Open an issue on GitHub for questions/bugs or contact me at kai.mason@ucalgary.ca for any questions about the work.
