# Privacy-Preserving Multi-Center Clinical Analytics with Threshold Homomorphic Encryption

This repository provides the implementation accompanying the paper *Privacy-Preserving Multi-Center Clinical Analytics with Threshold Homomorphic Encryption: A Pediatric Craniospinal Irradiation Use Case*.

The library implements three standard clinical statistical tests under the CKKS homomorphic encryption scheme, supporting both single-party and threshold (multi-party) settings:

- **Pearson correlation coefficient** -- for continuous dose-toxicity associations
- **Wilcoxon rank-sum test** -- for dose comparisons across toxicity groups
- **Chi-squared test** -- for categorical associations (e.g., radiation technique vs. lymphopenia)

Our code is built on top of the [OpenFHE](https://github.com/openfhe-development/openfhe-development) library. The ranking, order statistics, sorting, and median primitives are based on the [openfhe-statistics](https://github.com/FedericoMazzone/openfhe-statistics) library by Mazzone et al.


## Paper and Citation

If you use this work, please cite:

```bibtex
@inproceedings{mazzone2025efficient,
  title={Efficient Ranking, Order Statistics, and Sorting under CKKS},
  author={Mazzone, Federico and Everts, Maarten and Hahn, Florian and Peter, Andreas},
  booktitle={34th USENIX Security Symposium (USENIX Security '25)},
  year={2025},
  address={Seattle, WA},
  publisher={USENIX Association},
  month={aug}
}
```

and:

```bibtex
@article{turcas2025fhe,
  title={Privacy-Preserving Multi-Center Clinical Analytics with Threshold Homomorphic Encryption: A Pediatric Craniospinal Irradiation Use Case},
  author={{\c{T}}urca{\c{s}}, George C{\u{a}}t{\u{a}}lin and Gugulea, George and Lupa{\c{s}}cu, Cristian and Togan, Mihai and {\c{T}}urca{\c{s}}, Andrada Crina},
  year={2025}
}
```


## Contact

For bug reports or inquiries, please open an [issue](https://github.com/certfhe/fhe_medical_insights/issues).


## Prerequisites

### Compiler and Build Tools

```bash
sudo apt-get install build-essential cmake
```

### OpenFHE Library

This project requires [OpenFHE](https://github.com/openfhe-development/openfhe-development) **version 1.1.2**. Follow the steps below to install it:

1. **Clone the OpenFHE repository**

   ```bash
   git clone --branch v1.1.2 https://github.com/openfhe-development/openfhe-development.git
   cd openfhe-development
   ```

2. **Build and install**

   ```bash
   mkdir build && cd build
   cmake ..
   make -j$(nproc)
   sudo make install
   cd ../..
   ```

   If you do not have sudo access, specify a local installation prefix:

   ```bash
   cmake -DCMAKE_INSTALL_PREFIX=$HOME/openfhe ..
   make -j$(nproc)
   make install
   ```

   Then add the following to your shell configuration (`~/.bashrc` or `~/.zshrc`):

   ```bash
   export LD_LIBRARY_PATH=$HOME/openfhe/lib:$LD_LIBRARY_PATH
   export CMAKE_PREFIX_PATH=$HOME/openfhe:$CMAKE_PREFIX_PATH
   ```

For additional information, refer to the [OpenFHE documentation](https://openfhe-development.readthedocs.io/en/latest/).


## Building the Project

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
cd ..
```

If OpenFHE is installed in a non-standard location, specify it with:

```bash
cmake -DOpenFHE_DIR=/path/to/openfhe/build ..
```

After compilation, the following executables are available in `build/`:

| Executable | Description |
|---|---|
| `statistics-test` | Single-party Pearson, Wilcoxon, and Chi-squared tests |
| `mpc-statistics` | Multi-party (threshold) Pearson, Wilcoxon, and Chi-squared tests |
| `ranking` | Ranking benchmarks |
| `minimum` | Minimum-finding benchmarks |
| `median` | Median computation benchmarks |
| `sorting` | Sorting benchmarks |

### Running

```bash
./build/statistics-test
./build/mpc-statistics <num_values> <num_parties> [num_runs]
```


## Benchmarking

Each benchmark executable accepts the following arguments:

```bash
./build/ranking <vector_length> [<tie_correction>] [<single_thread>]
./build/minimum <vector_length> [<single_thread>]
./build/median  <vector_length> [<single_thread>]
./build/sorting <vector_length> [<single_thread>]
```

- **vector_length**: Length of the input vector (a power of two).
- **tie_correction** (ranking only): `1` to enable, `0` to disable (default: `0`).
- **single_thread**: `1` for single-threaded mode, `0` for multi-threaded (default: `0`).

### Automated Benchmarking

```bash
sh ./benchmark.sh <algorithm>
```

Where `<algorithm>` is one of: `ranking`, `ranking-tie`, `minimum`, `median`, `sorting`.

Results are saved to `benchmark.out`; logs are stored in `logs/`.


## Suggested Parameters

The parameter settings used for approximating the comparison and indicator functions can be found in the `test-*.cpp` files. We suggest using the **f, g approximation method** with the following composition degrees.

**Comparison function:**

| Bit Precision | df | dg |
|---|---|---|
| 1 | 2 | 1 |
| 2 | 2 | 1 |
| 3 | 2 | 2 |
| 4 | 2 | 2 |
| 5 | 2 | 3 |
| 6 | 2 | 3 |
| 7 | 2 | 4 |
| 8 | 2 | 4 |
| 9 | 2 | 5 |
| 10 | 2 | 5 |
| 11 | 2 | 5 |
| 12 | 2 | 6 |
| 13 | 2 | 6 |
| 14 | 2 | 7 |
| 15 | 2 | 7 |
| 16 | 2 | 8 |

**Indicator function:**
- df = 2
- dg = (log2(vectorLength) + 1) / 2


## Docker

A Dockerfile is provided for reproducibility:

```bash
docker build -t fhe-medical-insights .
docker run -it fhe-medical-insights
```

For optimal performance, we recommend using the standard (non-Docker) installation, as Docker can introduce runtime overhead for small vectors.


## Known Issues

**Missing shared library:**

If you encounter `error while loading shared libraries: libOPENFHEbinfhe.so.1`, the OpenFHE shared libraries are not in the linker search path. Fix with:

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

To make this permanent, add the line to `~/.bashrc` and run `source ~/.bashrc`.


## License

This project is licensed under the BSD 2-Clause License. See [LICENSE](LICENSE) for details.
