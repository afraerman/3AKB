[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C++-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

# 3AKB 🌍 — Three Algorithms for Gravitational Acceleration Benchmarking

**3AKB** is a lightweight, header-only C++20 console application for computing Earth's gravitational acceleration using high-degree spherical harmonic models (EGM96 and EGM2008).

It provides a unified framework to evaluate and benchmark three classical recursive algorithms — Belikov, Cunningham, and Holmes (Stokes) — under identical conditions, including a **multi-threaded parallel implementation** ⚡ of the Holmes algorithm described in our forthcoming paper.

**Key features**
- 🚀 Interactive console menu for quick experiments
- 📍 Single-point evaluation
- 🔄 Random-point and sequential orbit propagation comparison modes
- 📊 Comprehensive multithreaded benchmarking with speedup/efficiency metrics
- 📄 CSV export for post-processing and research workflows
- ✅ No external dependencies

## Table of Contents
- [Quickstart 🚀](#quickstart)
- [Installation 🛠️](#installation)
- [Project structure 📂](#project-structure)
- [Usage 🎮](#usage)
- [Input parameters ⚙️](#input-parameters)
- [Output format 📈](#output-format)
- [Algorithms / API reference 🔧](#algorithms--api-reference)
- [Modes 🧪](#modes)
- [Authors 👥](#authors)
- [Citation 📚](#citation)
- [License 📄](#license)

## Quickstart 🚀

1. Clone the repository:
   ```bash
   git clone https://github.com/timurezy/3AKB.git
   cd 3AKB
   ```
2. Open the Visual Studio solution `3AKB.sln`.
3. Build in **Release** configuration (x64 recommended).
4. Run the executable.

You will be greeted by an interactive menu:
```
1. RUN INDIVIDUAL ALGORITHM
2. ALGORITHM COMPARISON MODE
3. ALGORITHM BENCHMARKING MODE
4. CHANGE MAX HARMONICS
5. CHANGE INPUT COORDINATES
6. SELECT GRAVITY MODEL
7. IMPORT HARMONICS
8. NUMBER OF THREADS
0. EXIT
```

**Example**: select EGM2008, set NMAX = 2000, 12 threads → run benchmarking → get detailed CSV with speedup up to ~8.4× ⚡

## Installation 🛠️

### Requirements
- Windows 10/11
- Visual Studio 2019 or 2022 (with C++20 support)
- MSVC toolset

### Build
Open `3AKB.sln` → Build → Build Solution.  
No external libraries or package managers are required.

## Project structure 📂

**Core files**
- `3AKB.cpp` – program entry point and interactive menu
- `alltypes.h` – shared constants, gravity model selection, harmonic storage
- `SingleAlgorithmExecution.h` – single-point evaluation
- `AlgorithmComparison.h` – comparison mode + CSV export
- `AlgorithmBenchmarking.h` – multithreaded benchmarking
- `Simulate.h` – coordinate utilities and random point generation

**Algorithm implementations**
- `Belikov.h` – Belikov method
- `Cunningham.h` – Cunningham method
- `Albert.h` – single- and multi-threaded Holmes (Stokes) implementation

**Data**
- `EGM96.dat`, `EGM2008.dat` – harmonic coefficient files

## Usage 🎮

The program starts with an interactive console menu (see Quickstart).  
Note:
- `NMAX = 0` → central-body (Keplerian) acceleration only
- `NMAX > 0` → harmonic coefficients must be imported first

## Input parameters ⚙️

- **Gravity model**: EGM96 or EGM2008
- **Maximum degree/order (NMAX)**: 0 (Keplerian) to full model resolution
- **Point coordinates**: radius (m), geocentric latitude (°), longitude (°)
- **Number of threads**: 1–maximum supported by hardware (affects parallel Holmes only)
- **Number of runs**: used in comparison and benchmarking modes

## Output format 📈

All results are exported as CSV files in the executable directory.

**Random comparison** – `results_RANDOM.csv`  
Header: `Algorithm,Run,Time (ms),R,Latitude,Longitude,ax,ay,az`  
Includes summary statistics (total time, speedup, efficiency).

**Sequential propagation** – `results_SEQUENTIAL.csv`  
Similar format + final propagated position.

**Benchmarking** – multiple CSV files with detailed timing, speedup, and efficiency metrics.

## Algorithms / API reference 🔧

### Belikov
```cpp
void gravityBelikov(double r, double lat, double lon, int nmax,
                    std::array<double,3>& result);
```

### Cunningham
```cpp
void gravityCunningham(double r, double lat, double lon, int nmax,
                       std::array<double,3>& result);
```

### Holmes / Stokes (parallel-capable) ⚡
```cpp
namespace uniorb {
class gravity_stokes {
public:
    void use_concurrency(int thread_count);
    void get_acceleration(double r, double lat, double lon,
                          int nmax, std::array<double,3>& result);
};
}
```

## Modes 🧪

- **Individual algorithm**: single-point evaluation of selected method
- **Comparison**: random points or simple sequential propagation, with CSV export
- **Benchmarking**: systematic performance tests across degrees and thread counts

## Authors 👥

A. V. Fraerman¹, T. S. Pavlov¹, A. R. Shaykhutdinov¹, P. R. Zapevalin²  
¹ Astro Space Center, Lebedev Physical Institute, Russian Academy of Sciences  
² Sternberg Astronomical Institute, Lomonosov Moscow State University  

Corresponding author: fraerman@asc.rssi.ru

## Citation 📚

If you use this code in your research, please cite our paper "Multi-threaded CPU implementation of the Holmes algorithm for pointwise high-degree gravitational potential and acceleration evaluation":

```bibtex
@article{FRAERMAN2026,
title = {Multi-threaded CPU implementation of the Holmes algorithm for pointwise high-degree gravitational potential and acceleration evaluation},
journal = {Advances in Space Research},
volume = {78},
number = {5},
pages = {5267-5278},
year = {2026},
issn = {0273-1177},
doi = {https://doi.org/10.1016/j.asr.2026.06.089},
url = {https://www.sciencedirect.com/science/article/pii/S0273117726009154},
author = {A.V. Fraerman and T.S. Pavlov and A.R. Shaykhutdinov and P.R. Zapevalin},
keywords = {85–04},
abstract = {Accurate and efficient evaluation of Earth’s gravitational potential using high-degree spherical harmonic expansions is essential for precise long-term orbit propagation in celestial mechanics. Widely used recursive algorithms, such as those developed by Cunningham, Belikov, and Holmes, become computationally demanding at ultra-high degrees required by modern global gravity models. This paper introduces a multi-threaded parallel implementation of Holmes’ algorithm, designed for evaluations on standard multi-core CPUs. The method was extensively tested on various processor architectures and across a range of maximum harmonic degrees up to 2000. By distributing the outer order-wise summation across independent threads, the approach gives speedups of up to 8.4 times. Numerical tests confirm that the parallel version fully preserves the accuracy and stability of the serial Holmes algorithm. All algorithms discussed, including the parallel implementation, are freely available on GitHub.}
}
```

## License 📄

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
