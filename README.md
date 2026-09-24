# Privacy-Preserving Multi-Center Clinical Analytics with Threshold Homomorphic Encryption

This repository provides the implementation accompanying the paper *Privacy-Preserving Multi-Center Clinical Analytics with Threshold Homomorphic Encryption: A Pediatric Craniospinal Irradiation Use Case*.

The library implements three standard clinical statistical tests under the CKKS homomorphic encryption scheme, supporting both single-party and threshold (multi-party) settings:

- **Pearson correlation coefficient** -- for continuous dose-toxicity associations
- **Wilcoxon rank-sum test** -- for dose comparisons across toxicity groups
- **Chi-squared test** -- for categorical associations (e.g., radiation technique vs. lymphopenia)

Our code is built on top of the [OpenFHE](https://github.com/openfheorg/openfhe-development) library. The ranking, order statistics, sorting, and median primitives are based on the [openfhe-statistics](https://github.com/FedericoMazzone/openfhe-statistics) library by Mazzone et al.


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


## Changes for the Revision

- The masks of the client-aided step are now random: a new mask `u = exp(Z)`, with `Z` standard normal truncated to `|Z| <= 4`, is drawn for every query (`src/mask.h`).
- Wilcoxon and chi-squared now also end with a masked release, like Pearson.
- The two programs can read a cohort file (`--csv`) and pad it to a power of two. `mpc-statistics` splits it over the emulated centers, and each center encrypts its own rows.
- New command line options for the cohort size, the number of runs and the CKKS parameters.
- Each test has its own CKKS context (see below), the same in both programs.


## Contact

For bug reports or inquiries, please open an [issue](https://github.com/certfhe/fhe_medical_insights/issues).


## Prerequisites

### Compiler and Build Tools

```bash
sudo apt-get install build-essential cmake
```

On macOS, install `cmake` and `libomp` with Homebrew (`brew install cmake libomp`).

### OpenFHE Library

This project requires [OpenFHE](https://github.com/openfheorg/openfhe-development) **version 1.1.2**. Follow the steps below to install it:

1. **Clone the OpenFHE repository**

   ```bash
   git clone --branch v1.1.2 https://github.com/openfheorg/openfhe-development.git
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
# threshold version, 3 parties, synthetic cohort of 64 patients, 50 runs
./build/mpc-statistics --parties 3 --n 64 --runs 50

# threshold version on a cohort file, 2 parties
./build/mpc-statistics --parties 2 --csv data/example_cohort.csv --runs 1

# single-party version
./build/statistics-test --n 64 --runs 50
./build/statistics-test --csv data/example_cohort.csv --runs 1
```

Each program prints, for every analysis, the plaintext value, the value obtained from the encrypted computation, the deviation and the average time of one encrypted evaluation (including the masked release and the decryption). The last block, `SUMMARY`, has one line per analysis.

Main options (the same in both programs, except `--parties`):

- **--parties 2|3**: Number of key-share holders (`mpc-statistics` only, default: `3`).
- **--n N**: Size of the synthetic cohort (default: `512`). Ignored with `--csv`.
- **--runs R**: Number of timed runs (default: `50`).
- **--csv FILE**: Run the five analyses of the paper on a cohort file (format below).
- **--no-wilcoxon**: Skip the Wilcoxon test. It is skipped anyway for more than 64 values, since the ranking needs n^2 slots.
- **--verbose**: Print the CKKS parameters, the key sizes and, for every run, the masked values and the statistic computed from them.

The CKKS parameters can be changed with `--depth`, `--intbits`, `--decbits` (Pearson), `--chi2-depth`, `--chi2-intbits`, `--chi2-decbits` (chi-squared) and `--compare-depth`, `--rank-intbits`, `--rank-decbits` (ranking; the depth of the ranking context is set with `--wilcoxon-depth`, in `mpc-statistics` only). The defaults are the ones used in the paper:

| Test | Depth | Bits (integer + fractional) | Ring dimension |
|---|---|---|---|
| Pearson | 8 | 1+35 | 16384 |
| Chi-squared | 5 | 10+50 | 16384 |
| Wilcoxon (ranking) | 14 (12 in `statistics-test`) | 1+40 | 32768 |

In `mpc-statistics` the rank vector is refreshed with interactive bootstrapping before the final step. `statistics-test --bootstrap` uses ordinary CKKS bootstrapping instead of the deeper context (not used for the paper). The synthetic cohort is generated from a fixed seed, but `std::normal_distribution` is implemented differently in libstdc++ and libc++, so the synthetic values on Linux are not the same as on macOS (the tables of the paper were made on macOS).

### Masked Release

All three tests end with a short client-aided step. The server multiplies a few aggregates by fresh random masks, only these masked values are (threshold) decrypted, and the analyst computes the statistic from them. The masks cancel:

- **Pearson**: `(u1 XdotY/n, u2^2 LenX/n, u3^2 LenY/n)` with `u3 = u1/u2`, then `r = m1 / (sqrt(m2) sqrt(m3))`.
- **Wilcoxon**: `(u D1/n^2, u^2 D2/n^4)` with `D1 = S1 - n1(n+1)/2` and `D2 = n0 n1 (n+1)/12`, then `z = m1 / sqrt(m2)`.
- **Chi-squared**: `(v (ad-bc)^2/N^4, v r0 r1 c0 c1/N^4)`, whose ratio is `phi^2 = chi2/N`, then `chi2 = N m1/m2`.

The divisions by powers of n keep the released values small; without them the CKKS decoder can fail for larger cohorts.

### Cohort File

`--csv` reads a comma-separated file with a header line and one row per patient (an empty field is a missing value):

```
patient_index,center,pelvic_dmean_gy,cervical_dmean_gy,alc_week1,anc_week1,lymphopenia_grade_max,neutropenia_grade_max,technique_vmat_ht
```

The five analyses of Table 1 of the paper are computed on complete cases: two Pearson correlations (pelvic dose vs. lymphocytes, cervical dose vs. neutrophils, week 1), two Wilcoxon tests (pelvic dose by grade >= 3 lymphopenia, cervical dose by grade >= 3 neutropenia) and one chi-squared test (technique vs. grade 4 lymphopenia). `mpc-statistics` assigns the rows to the emulated centers in turn (row i of the file goes to center i mod parties). The `center` column has to be there, but it is not used.

`data/example_cohort.csv` is a **synthetic** file with the same shape as the clinical cohort (51 rows, with missing values). There is **no patient data** in this repository. The runs on the clinical cohort were done inside the originating institution with this code, and only the aggregate results were taken out.


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

**Memory usage of `mpc-statistics`:**

The threshold program emulates all parties in one process, so it keeps the key shares and evaluation keys of every party in memory at the same time. With the Wilcoxon test (ranking context, ring dimension 32768) the peak memory for a vector length of 64 was about 6.6 GB with 2 parties and 9.1 GB with 3 parties on macOS. On machines with less memory, use `--no-wilcoxon`, or `statistics-test` for the Wilcoxon test (about 1.4 GB).


## License

This project is licensed under the BSD 2-Clause License. See [LICENSE](LICENSE) for details.
