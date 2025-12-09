# AI Coding Agent Instructions for `qtlhot`

## Overview
The `qtlhot` package is an R-based library for statistical inference of QTL (Quantitative Trait Loci) hotspots and causal modeling in system genetics. It includes functions for co-mapping trait hotspots and causal model inference. The package is based on methodologies described in:

- Chaibub Neto et al. (2012): Quantile-based permutation thresholds for QTL hotspots. [Genetics, 191:1355-1365](https://doi.org/10.1534/genetics.112.139451).
- Chaibub Neto et al. (2013): Modeling causality for pairs of phenotypes in system genetics. [Genetics, 193:1003-1013](https://doi.org/10.1534/genetics.112.147124).

## Codebase Structure
- **R/**: Contains the core R scripts implementing the package's functionality.
- **man/**: Documentation files for the R functions.
- **data/**: Example datasets used in the package.
- **vignettes/**: Tutorials and examples in R Markdown format.
- **inst/notes/**: Legacy and exploratory scripts, including deprecated code.
- **parallel/**: Scripts for parallelized workflows.

## Key Patterns and Conventions
1. **Function Documentation**:
   - All functions are documented in the `man/` directory.
   - Follow the Roxygen2 format for inline documentation.

2. **Data Handling**:
   - Example datasets are stored in `data/` as `.RData` files.
   - Use `load()` to access these datasets in scripts.

3. **Testing and Validation**:
   - Statistical methods often rely on permutation tests (e.g., `hotperm.R`).
   - Ensure reproducibility by setting random seeds where applicable.

4. **Parallelization**:
   - Scripts in `parallel/` demonstrate how to distribute computations across multiple cores.
   - Refer to `runparallel.sh` for an example workflow.

5. **Vignettes**:
   - Tutorials in `vignettes/` provide usage examples.
   - Use `knitr` and `rmarkdown` for rendering.

## Developer Workflows

Ensure that Pandoc is installed to rebuild vignettes. On macOS, install it using Homebrew:
```zsh
brew install pandoc
```

### Building the Package

Run the following command in the R console to build the package:
```R
R CMD build .
```

### Installing the Package

Install the package locally using:

```R
R CMD INSTALL qtlhot_1.2.3.tar.gz
```

### Testing
- Use the example datasets in `data/` to validate functions.
- Check the vignettes for end-to-end examples.

### Debugging
- Use `browser()` within R scripts to set breakpoints.
- Leverage `traceback()` to debug errors.

### Modifying Code
- Do not alter NAMESPACE, DESCRIPTION or *.Rd files manually; use Roxygen2 comments and `devtools::document()` instead.

### Running Package Checks

After a round of development, follow these steps to run package checks:

- Commit and push all changes to the repository.
- Update documentation files using `devtools::document()`.
- Build the repository package with `devtools::build()`.
- Install the built package using `devtools::install()`.
- Run package checks with `devtools::check(document = FALSE, args = c('--as-cran'))`.

```zsh
R -e "devtools::document()"
R -e "devtools::build()"
R -e "devtools::install()"
R -e "devtools::check(document = FALSE, args = c('--as-cran'))"
```

The `devtools::check()` function is used to ensure the package meets CRAN standards. The options used here are:
- **`document = FALSE`**: Skips re-documenting the package.
- **`args = c('--as-cran')`**: Runs checks as if submitting to CRAN.

The check command will:
- Build the package.
- Check for errors, warnings, and notes.
- Validate examples, vignettes, and documentation.

### Fixing Issues from Package Checks
- Address any errors, warnings, or notes reported by `devtools::check()`.

## External Dependencies
- R packages: `stats`, `qtl`, `mnormt`, `utils`, `corpcor`, `broman`.
- Suggested packages: `knitr`, `rmarkdown` for vignettes.

## Notes
- Legacy code and deprecated scripts are stored in `inst/notes/`.
- The `DESCRIPTION` file provides metadata and dependency information.
- For questions, contact the maintainer: Brian S. Yandell (<brian.yandell@wisc.edu>).

## Examples
### Running a Permutation Test
```R
load("data/hotperm1.RData")
result <- hotperm(data)
plot(result)
```

### Generating a Vignette
```R
library(knitr)
knit("vignettes/qtlhot.Rmd")
```