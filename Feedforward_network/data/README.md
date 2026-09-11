# Data Folder

Dataset files are not distributed with this repository. Download them from the
upstream sources and prepare the required MAT files locally using
`../docs/DATASETS.md`, which records the sources, licences, encodings, variable
names, and array formats. Files placed here are ignored by Git. The repository's
MIT licence does not replace the upstream dataset licences documented there and
in `../../THIRD_PARTY_NOTICES.md`.

Required dataset files for a full local run are:

- `abalone_dataset.mat`
- `breast_cancer_dataset.mat`
- `car_dataset.mat`
- `mushroom_dataset.mat`
- `toyota_dataset.mat`
- `mnist.mat`
- `mnist_fashion.mat`
- `kmnist.mat`
- `afro_mnist_ethiopic.mat`
- `afro_mnist_nko.mat`
- `afro_mnist_osmanya.mat`
- `afro_mnist_vai.mat`

`fisheriris` is loaded directly from MATLAB and does not belong in this folder.
Pong supervised data and
motor-control episodes are generated inside their respective training scripts.
Motor-control training writes
`outputs/training_data/motor_control/episodes_seed_###.mat` as a generated
artefact; it is not a required input.
