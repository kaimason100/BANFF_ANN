# Dataset Preparation

This document describes every dataset file expected by the active scripts. Put
all prepared dataset files in the repository `data` folder.

The code does not download datasets automatically. This is intentional: it keeps
the training scripts deterministic and avoids silently changing source data.

## Required files

The active local package expects these files when running all tasks:

```text
data/abalone_dataset.mat
data/breast_cancer_dataset.mat
data/car_dataset.mat
data/mushroom_dataset.mat
data/toyota_dataset.mat
data/mnist.mat
data/mnist_fashion.mat
data/kmnist.mat
data/afro_mnist_ethiopic.mat
data/afro_mnist_nko.mat
data/afro_mnist_osmanya.mat
data/afro_mnist_vai.mat
data/LQR_SUPERVISED_EPISODES_20251002_172744.mat
examples/dynamical_systems/dynamics_list.xlsx
```

`fisheriris` is loaded directly from MATLAB and does not require a file in
`data`.

Pong supervised samples are generated inside the Pong training script and do not
require an external dataset.

## Tabular classification datasets

### Breast cancer

Expected file:

```text
data/breast_cancer_dataset.mat
```

Expected variable:

```matlab
data
```

Supported formats:

- table: column 2 is the class label, columns 3:end are numeric features
- cell array: column 2 is the class label, columns 3:end are numeric features
- numeric matrix: column 2 is the class label, columns 3:end are numeric
  features
- struct with fields `features` and `labels`, where rows are samples

Common source:

- UCI Breast Cancer Wisconsin Diagnostic:
  https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic

If using the UCI diagnostic CSV, remove or ignore the ID column, keep diagnosis
as the label, and put numeric features after the label. Save as:

```matlab
save("data/breast_cancer_dataset.mat", "data")
```

### Car quality

Expected file:

```text
data/car_dataset.mat
```

Expected variable:

```matlab
data
```

Expected shape:

```matlab
N x 7 numeric matrix
```

Columns 1:6 are features. Column 7 is the class label. The loader converts
column 7 to `categorical`.

Common source:

- UCI Car Evaluation:
  https://archive.ics.uci.edu/dataset/19/car+evaluation

The raw UCI file is categorical text. Convert each categorical column to a
consistent numeric code before saving. For example, after reading the CSV into a
table `T`, convert each column with `grp2idx(categorical(...))`, concatenate the
columns, and save the numeric matrix as `data`.

### Mushroom

Expected file:

```text
data/mushroom_dataset.mat
```

Expected variable:

```matlab
data
```

Expected shape:

```matlab
N x M numeric matrix
```

Column 1 is the class label. Columns 2:end are features. The loader converts
column 1 to `categorical`.

Common source:

- UCI Mushroom:
  https://archive.ics.uci.edu/dataset/73/mushroom

The raw UCI file is categorical text. Convert each categorical column to a
consistent numeric code before saving.

## Regression datasets

### Abalone

Expected file:

```text
data/abalone_dataset.mat
```

Expected variable:

```matlab
data
```

Expected shape:

```matlab
N x M numeric matrix
```

Columns 1:end-1 are features. The final column is the target.

Common source:

- UCI Abalone:
  https://archive.ics.uci.edu/dataset/1/abalone

The original UCI file has a categorical sex column followed by numeric
measurements and rings. Convert sex to a numeric code or one-hot columns before
saving. If reproducing the current loader exactly, the target must be the final
column.

### Toyota

Expected file:

```text
data/toyota_dataset.mat
```

Expected variable:

```matlab
data
```

Expected shape:

```matlab
N x M numeric matrix
```

The current loader uses:

- target: column 3
- features: columns 1, 2, and 4:end

This means column 3 must be the regression target for this package. If preparing
data from a Toyota Corolla CSV whose target is named `Price`, arrange the saved
numeric matrix so that `Price` is column 3, and put all selected numeric or
encoded predictor columns in the remaining columns.

Common public versions of this dataset are the Toyota Corolla used-car price
data, often described as 1436 records of used Toyota Corollas sold in the
Netherlands with attributes such as price, age, kilometers, horsepower, doors,
tax, and weight.

## MNIST-family image datasets

Expected files:

```text
data/mnist.mat
data/mnist_fashion.mat
data/kmnist.mat
data/afro_mnist_ethiopic.mat
data/afro_mnist_nko.mat
data/afro_mnist_osmanya.mat
data/afro_mnist_vai.mat
```

Each file must contain:

```matlab
training
test
```

The loader expects `training` and `test` to contain image data and labels. The
supported internal field names are handled by the script, but the safest format
is:

```matlab
training.images   % 28 x 28 x 1 x N, or 28 x 28 x N, or 784 x N
training.labels   % N x 1 labels, values 0:9 or categorical labels
test.images       % same image layout for held-out test data
test.labels       % test labels
```

Images should be grayscale 28 x 28. Labels should identify ten classes. The
training script creates its own train/validation split from `training`; it uses
the supplied `test` split only for final testing.

Common sources:

- MNIST: http://yann.lecun.com/exdb/mnist/
- Fashion-MNIST: https://github.com/zalandoresearch/fashion-mnist
- Kuzushiji-MNIST: https://github.com/rois-codh/kmnist
- Afro-MNIST: https://github.com/Daniel-Wu/AfroMNIST

Save each converted dataset as a MAT file with the expected filename. Example
shape check:

```matlab
load("data/mnist.mat")
size(training.images)
size(training.labels)
size(test.images)
size(test.labels)
```

## Motor-control data

Expected file:

```text
data/LQR_SUPERVISED_EPISODES_20251002_172744.mat
```

The active training script generates its supervised episodes internally and also
uses the LQR setup functions embedded in the live script. The MAT file is kept
for compatibility with existing plotting and saved-data workflows. Do not
rename it unless you also update the scripts that reference it.

## Dynamical-system definitions

Expected file:

```text
examples/dynamical_systems/dynamics_list.xlsx
```

The active DS trainer reads this spreadsheet. Required columns are accessed by
position:

- column 1: system name
- column 2: state dimension
- columns 4:6: initial-condition entries

The active script uses task rows:

```matlab
[5, 34, 39, 41, 47, 32, 8, 9, 25]
```

Those rows correspond to:

```text
Lorenz, MO0, MO5, MO7, MO13, Rikitake, SprottB, SprottC, SprottS
```

If editing `dynamics_list.xlsx`, keep the row ordering and columns stable unless
you intentionally update the task-row list in the training and test scripts.

## Minimal conversion examples

Convert a numeric CSV whose final column is the target:

```matlab
T = readtable("raw_file.csv");
D = table2array(T);
data = double(D);
save("data/abalone_dataset.mat", "data")
```

Convert categorical text columns to numeric codes:

```matlab
T = readtable("raw_categorical_file.csv", TextType="string");
data = zeros(height(T), width(T));
for c = 1:width(T)
    if isnumeric(T{:, c})
        data(:, c) = double(T{:, c});
    else
        data(:, c) = grp2idx(categorical(T{:, c}));
    end
end
save("data/car_dataset.mat", "data")
```

Create the recommended MNIST-family structure after loading arrays:

```matlab
training.images = XTrain;   % 28 x 28 x 1 x N
training.labels = YTrain;   % N x 1
test.images = XTest;
test.labels = YTest;
save("data/mnist.mat", "training", "test", "-v7.3")
```

## Data-leakage rules

Do not normalize or standardize the full dataset before saving it. Save raw
numeric features and labels. The training scripts fit normalization statistics
from the training split only, save those statistics with each network, and the
tests assert that the saved statistics match the training split.

Do not merge validation or test examples into the training split. The package
creates fixed train/validation/test splits internally for tabular tasks and
uses the official held-out `test` split for MNIST-family tasks.

