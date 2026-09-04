# Data Folder

This folder contains the exact prepared MAT datasets used by the released local
training and test scripts. They are supplied with the repository, so the normal
workflow does not require a separate dataset download. Verify the local files
with `SHA256SUMS`. See
`../docs/DATASETS.md` for sources, licences, exact encodings, expected variable
names, and array formats. The repository's MIT licence does not replace the
third-party dataset licences documented there and in
`../../THIRD_PARTY_NOTICES.md`.

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

`fisheriris` is loaded directly from MATLAB. Pong supervised data and
motor-control episodes are generated inside their respective training scripts.
Motor-control training writes
`LQR_SUPERVISED_EPISODES_20251002_172744.mat` as a generated compatibility
artefact; it is not a required input and is not included in the public dataset
bundle.
