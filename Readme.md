## Abstract

We present a user-friendly, standalone C/C++ program for generating molecular potential energy
curves (PECs) using the Rydberg--Klein--Rees (RKR) method and experimentally determined
spectroscopic constants. The program uses the open-source GNU Scientific Library (GSL) for
numerical computations and provides a graphical user interface (GUI) to enter constants and
calculation parameters. It also visualizes the resulting potential energy curve in the same interface
for quick verification and exports it to a plain-text format compatible with the LEVEL 8.0 program,
enabling downstream calculations such as Franck-Condon factors for comparing the rovibrational
peaks with the experimentally obtained spectra. The program is benchmarked against a proof-of-concept Matlab (MathWorks, Natick, MA) script developed by our group as well as the RKR script in the pre-existing
Optimizer Matlab package developed by Sovkov (https://sourceforge.net/projects/optimizer-sovkov/). It shows significant improvement in compute-time while maintaining
the same double-precision (64-bit floating-point) numerical accuracy. Most importantly, the program
is free-to-use under the GNU General Public License and is purely standalone. Therefore, it is
accessible to universities and colleges without relying on proprietary software or requiring a
complex installation and external dependencies.


## Introduction
This is a C++ project for Rydberg-Klein-Reese (RKR) Procedure for point-wise diatomic molecular potential energy curves (PECs). The standalone executable shipped with this project **RKR_v1.exe** (compiled using msys/mingw64 build system, further information and installation guidelines can be found here: https://www.msys2.org/) offers fast, easy-to-use, and accessible (completely free) gui for inputing spectroscopic parameters/constants obtained from experimental data (laser-induced fluorescence spectra), computational backend generating the RKR PECs, plotting for immediate visualization, and saving the PEC data for further spectroscopy and quantum mechanical applications highlighted in the abstract above. 

---

## License

This project is licensed under the **GNU General Public License v3.0 or later**.  
See the `LICENSE` file for details.

### Third-Party Dependencies

| Library | Location | License |
|-----------|----------|----------|
| Eigen | `lib/eigen-master/` | Mozilla Public License v2.0 |
| GNU Scientific Library (GSL) | `lib/gsl/` | GNU General Public License v3.0 or later |

---

## Requirements
msys/mingw64 build system or similar C/C++ build environment.
https://www.msys2.org/


### GUI (Windows)

Download and install **wxWidgets**:

https://wxwidgets.org/downloads/

---

## Libraries

### GNU Scientific Library (GSL)

C++ numerical library providing scientific computing routines.

The source code is included in this project under:

```
lib/gsl/
```

Alternatively, clone the official repository:

```bash
git clone git://git.savannah.gnu.org/gsl.git
```

---

### Eigen

A C++ template library for linear algebra, including matrices, vectors, numerical solvers, and related algorithms.

The source code is included in:

```
lib/eigen-master/
```

Alternatively, clone the repository:

```bash
git clone https://github.com/PX4/eigen.git
```

---

## Project Structure

```
.
├── src/                 # Source code
├── lib/
│   ├── gsl/             # GNU Scientific Library
│   └── eigen-master/    # Eigen linear algebra library
├── molecular_data/      # Molecular constants and input datasets
├── LICENSE
└── README.md
```

---

## Directory Description

- **src/**  
  Contains the C++ source files.

- **lib/gsl/**  
  GNU Scientific Library source code and build files.

- **lib/eigen-master/**  
  Eigen header-only linear algebra library.

- **molecular_data/**  
  Molecular constants and spectroscopy data files used by the program.

---

## Building

Before building, ensure:

1. **wxWidgets** is installed (for GUI support).
2. **GSL** is available (already included in `lib/gsl/`).
3. **Eigen** is needd to be git cloned (include in `lib/eigen-master/`).

Then compile the project using the provided Makefile:

```bash
make
```

---

## Authors

**Amar Dadel**, **Dr. Burcin Bayram**  
dadela@miamioh.edu, bayramsb@miamioh.edu \
Department of Physics, Miami University Oxford, Ohio, United States.

