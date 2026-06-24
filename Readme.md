# Project Name

A C++ project for molecular spectroscopy and scientific computation.

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
Department of Physics, Miami University Oxford, Ohio \ 
United States.

