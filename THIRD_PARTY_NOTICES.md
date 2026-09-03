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
