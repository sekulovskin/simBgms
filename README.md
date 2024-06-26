# simBgms: Simplifying Simulation Studies for Bayesian Markov Random Field Models 

## Description

The `simBgms` package in R provides a simplified approach to running simulation studies for Bayesian graphical modeling. Typically, researchers want to make decisions about sample size, number of observations, and the specific choices of prior distributions to use. With this package, users can easily simulate data for Markov random field models and then estimate these models using the `bgms` and `BDgraph` packages in R. The package streamlines the process of running simulation studies, eliminating the need for advanced programming skills. In addition, the package provides the ability to parallelize model estimation, which increases computational efficiency.

**Please note** that there has been an update to the options for simulating the data! Therefore, readers who come to this site as a result of reading Sekulovski et al. (2024) should refer to the updated documentation.

## Key Features

- Simulate data for Markov random field models.
- Estimate models using the `bgms` and `BDgraph` packages.
- Simplified approach for conducting simulation studies.
- Parallelize model estimation for improved computational efficiency.

## Installation

You can install the `simBgms` package from GitHub using the `devtools` package:

```R
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
devtools::install_github("sekulovskin/simBgms")
```

## Usage

```R
library(simBgms)
sim_bgm(level = "Discrete",
        variable_type = "ordinal",
        repetitions = 10, 
        no_observations = c(100, 500), 
        no_variables = c(5, 10),
        no_categories = c(1, 2), 
        ordinal_sim_type = "irt",
        induce_sparsity = TRUE,
        induce_sparsity = FALSE,
        density = c(0.1, 0.5, 0.9))
```

For detailed usage instructions and examples, please refer to the package documentation.

## Contributing

Contributions to `simBgms` are welcome! If you encounter any issues or have suggestions for improvement, please feel free to open an issue or submit a pull request on the GitHub repository.

## License

This package is licensed under the [MIT License](https://opensource.org/licenses/MIT). See the LICENSE file for more details.

## References

Sekulovski, N., Keetelaar, S., Haslbeck, J., & Marsman, M. (2024). Sensitivity analysis of prior distributions in Bayesian graphical modeling: Guiding informed prior choices for conditional independence testing. advances.in/psychology, 2, e92355. https://doi.org/10.56296/aip00016

