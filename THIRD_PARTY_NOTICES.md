# Third-Party Notices

This repository includes or interoperates with the third-party material listed
below. These notices apply only to the identified material; the remainder of
the repository is licensed under the MIT License in `LICENSE`.

## slanColor / slanCM

The feedforward publication plotting support includes a local `slanCM`
compatibility function. Its interface and fallback colour-map data are based
on colour maps distributed by the slanColor project:

- Project: slanColor (`slanCM`)
- Author: Zhaoxu Liu / slandarer
- Source: <https://github.com/slandarer/slanColor>
- MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormaps>
- Licence: BSD 3-Clause

The complete slanColor BSD 3-Clause notice is retained in
`LICENSES/BSD-3-Clause-slanColor.txt`.

The local compatibility function is not the complete upstream slanColor
package. If an independently installed copy of `slanCM` is available on the
MATLAB path, the plotting helper may call it; that installation remains subject
to its own accompanying licence and notices.

The upstream collection also acknowledges the original sources of individual
colour maps, including Matplotlib, scicomap, CMasher, Scientific Colour Maps,
cmocean, palettable and colorcet. See the upstream slanColor documentation for
the detailed provenance and citations for each map.

## Prepared Dataset Files

`Feedforward_network/data` contains prepared MAT-file copies of the datasets
listed below. Preparation comprised MATLAB-format conversion and, for tabular
categorical datasets, numeric category encoding. The files remain subject to
their upstream licences; the repository MIT licence does not supersede them.

### UCI datasets: CC BY 4.0

The following datasets are redistributed under the Creative Commons
Attribution 4.0 International licence:
<https://creativecommons.org/licenses/by/4.0/>.

- Abalone, UCI Machine Learning Repository,
  <https://doi.org/10.24432/C55C7W>
- Breast Cancer Wisconsin (Diagnostic), William Wolberg, Olvi Mangasarian,
  Nick Street and William Street, UCI Machine Learning Repository,
  <https://doi.org/10.24432/C5DW2B>
- Car Evaluation, Marko Bohanec, UCI Machine Learning Repository,
  <https://doi.org/10.24432/C5JP48>
- Mushroom, UCI Machine Learning Repository,
  <https://doi.org/10.24432/C5959T>

The prepared files are `abalone_dataset.mat`,
`breast_cancer_dataset.mat`, `car_dataset.mat`, and
`mushroom_dataset.mat`.

### Toyota used-car data: CC0 1.0

`toyota_dataset.mat` is derived from the Toyota subset of Aditya Desai's
100,000 UK Used Car dataset:
<https://www.kaggle.com/datasets/adityadesai13/used-car-dataset-ford-and-mercedes>.
The upstream dataset is marked CC0: Public Domain. See
<https://creativecommons.org/publicdomain/zero/1.0/>.

### MNIST: CC BY-SA 3.0

`mnist.mat` contains the MNIST database by Yann LeCun, Corinna Cortes and
Christopher J. C. Burges, obtained from
<https://yann.lecun.com/exdb/mnist/>. MNIST is made available under the
Creative Commons Attribution-ShareAlike 3.0 licence:
<https://creativecommons.org/licenses/by-sa/3.0/>.

### Fashion-MNIST: MIT

`mnist_fashion.mat` is derived from Fashion-MNIST, copyright 2017 Zalando SE:
<https://github.com/zalandoresearch/fashion-mnist>. The upstream MIT licence is
retained in `LICENSES/MIT-Fashion-MNIST.txt` and is available upstream at
<https://github.com/zalandoresearch/fashion-mnist/blob/master/LICENSE>.

### Kuzushiji-MNIST: CC BY-SA 4.0

`kmnist.mat` is derived from Kuzushiji-MNIST by ROIS-DS CODH:
<https://github.com/rois-codh/kmnist>. It is made available under the Creative
Commons Attribution-ShareAlike 4.0 International licence:
<https://creativecommons.org/licenses/by-sa/4.0/>.

### Afro-MNIST: MIT

`afro_mnist_ethiopic.mat`, `afro_mnist_nko.mat`,
`afro_mnist_osmanya.mat`, and `afro_mnist_vai.mat` are derived from
Afro-MNIST by Daniel J. Wu, Andrew C. Yang and Vinay U. Prabhu:
<https://github.com/Daniel-Wu/AfroMNIST>. The upstream MIT licence is available
at <https://github.com/Daniel-Wu/AfroMNIST/blob/master/LICENSE> and is retained
in `LICENSES/MIT-AfroMNIST.txt`.
