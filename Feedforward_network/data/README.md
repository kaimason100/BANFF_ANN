# Data Folder

Place local dataset MAT files here when running the training and test scripts.
The release zip intentionally does not include the large task datasets or saved
outputs. See `docs/DATASETS.md` for download sources, expected filenames,
variable names, and array/table formats.

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

`fisheriris` is loaded directly from MATLAB. Pong supervised data and the
motor-control LQR episodes are generated inside their respective training
scripts. Motor-control training writes
`LQR_SUPERVISED_EPISODES_20251002_172744.mat` as a generated compatibility
artifact; it is not a required input for a clean run.
