# Rydberg-Klein-Reese (RKR) Potential Energy Curve Generator

This is a C++ project for Rydberg-Klein-Reese (RKR) Procedure for point-wise diatomic molecular potential energy curves (PECs). The standalone executable shipped with this project **RKR_v1.exe** (compiled using msys/mingw64 build system, further information and installation guidelines can be found here: https://www.msys2.org/) offers fast, easy-to-use, and accessible (completely free) gui for inputing spectroscopic parameters/constants obtained from experimental data (laser-induced fluorescence spectra), computational backend generating the RKR PECs, plotting for immediate visualization, and saving the PEC data for further quantum mechanical applications, such as solving radial schrödinger equation

---

## License

This project is licensed under the **GNU General Public License v3.0**.  
See the `LICENSE` file for details.

### Third-Party Dependencies

| Library | Location | License |
|-----------|----------|----------|
| Eigen | `lib/eigen-master/` | Mozilla Public License v2.0 |
| GNU Scientific Library (GSL) | `lib/gsl/` | GNU General Public License v3.0 or later |

---

## Requirements

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

