# opticx

<div align=center>

![GitHub release (with filter)](https://img.shields.io/github/v/release/xatu-code/opticx)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/xatu-code/opticx/issues)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41524--024--01504--2-red.svg)](https://doi.org/10.1038/s41524-024-01504-2)
[![Documentation Status](https://readthedocs.org/projects/opticx/badge/?version=latest&style=flat)](https://opticx.readthedocs.io/en/latest/)

</div>

opticx computes the **linear and nonlinear optical response of crystals**, both at the
independent-particle level and including excitonic effects. Starting from a Wannier90 tight-binding
model (orthonormal, or with an overlap matrix) — and, optionally, exciton eigenstates obtained with
[Xatu](https://github.com/alejandrojuria/xatu) — it evaluates the first- and second-order optical
conductivities: absorbance, second-harmonic generation, the electro-optic (Pockels) effect, optical
rectification and the shift current. All of the second-order processes are branches of one general
two-frequency response σ(ω₁+ω₂; ω₁, ω₂), which can also be scanned as a full two-dimensional map.

<p align="center">
  <img src="hbn_opticx_example.png" width="90%" height="90%">
</p>

The theory behind the code, the conventions it follows and usage examples are described in the
[documentation](#documentation). The implementation follows
[Esteve-Paredes *et al.*, *npj Computational Materials* **11**, 13 (2025)](https://doi.org/10.1038/s41524-024-01504-2);
usage of the code requires citing the paper.

## Installation

opticx is written in Fortran and needs a Fortran compiler, OpenMP, and BLAS/LAPACK through OpenBLAS.

### Ubuntu 22.04 LTS native and WSL

Install the required libraries:
```
sudo apt-get install gfortran libopenblas-dev
```

Then build the binary (`bin/opticx`) with:
```
make
```

A run is launched by passing it an input file:
```
bin/opticx input.txt
```

### MacOS

The dependencies can be installed via `brew`. Use brew's `gcc` rather than the system compiler:
```
brew install gcc openblas
```

Then set the compiler and the library location in the Makefile:
```
FC     = gfortran-13
LIBS   = -L/opt/homebrew/opt/openblas/lib -lopenblas -fopenmp -lgfortran
```

### Intel MKL

The Makefile also carries an MKL variant of `LIBS`; uncomment it in place of the OpenBLAS one if you
would rather link against MKL.

## Tests

The test suite needs Python 3 with NumPy, and matplotlib if you want the diagnostic plots:
```
make test test_matrix                       # kernel equivalence
make run_test_shg_consistency               # SHG kernels on synthetic data
make run_test_shift_real run_test_shg_real  # real-data physics checks
make run_test_second_symmetry               # two-frequency symmetry checks
make check_sp_shift check_sp_shg            # against independent NumPy evaluations
make check_shift_covariant                  # degenerate-band (covariant) methods
make check_realtime_sign                    # absolute sign and normalisation vs a real-time simulation
make check_ex_rectification                 # excitonic rectification, injection current included
make check_out_of_plane                     # out-of-plane components and injection sign vs real time
make check_gauge_covariance                 # invariance under eigenvector gauge changes
make check_tb_hermiticity                   # model-file Hermiticity check and repair
make check_bandlist_guard                   # band-window report and unusual-Bandlist warning
make check_bands                            # band structure along a k-path (Kpath, Response = bands)
make check_ome_cache check_a4_basis_guard   # matrix-element cache and basis guards
```
Each prints `ALL TESTS PASSED` or `ALL CHECKS PASSED`.

## Documentation

The documentation is built with Sphinx and hosted on Read the Docs. It covers installation, the input
file reference, the output formats, the conventions the code follows (frequency axis, broadening,
units — worth reading before comparing numbers with a paper), the second-order theory, and the
validation suite.

The documentation lives in its own repository rather than in this one, so it can be edited and
rebuilt without touching the code. To build it locally, clone that repository and use a virtual
environment (it needs `sphinx-rtd-theme`, which most system Pythons do not have):
```
python3 -m venv .venv
.venv/bin/pip install -r docs/requirements.txt
.venv/bin/sphinx-build -b html docs docs/_build/html
```
and open `docs/_build/html/index.html`.

## License

opticx is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for the full text.
