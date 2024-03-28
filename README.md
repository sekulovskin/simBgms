# simBgms: Simplifying Simulation Studies for Bayesian Markov Random Field Models 

## Description

The `simBgms` package in R provides a simplified approach to conducting simulation studies for Bayesian graphical modeling. With this package, users can easily simulate data for Markov random field models and subsequently estimate these models using the `bgms` and `BDgraph` packages in R. It streamlines the process of performing simulation studies, eliminating the need for advanced programming skills. Additionally, the package offers the capability to parallelize model estimation, enhancing computational efficiency.

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
        induce_sparsity = FALSE,
        density = c(0.1, 0.5, 0.9), 
        irt = TRUE, 
        no_cores = 2)
```

For detailed usage instructions and examples, please refer to the package documentation.

## Contributing

Contributions to `simBgms` are welcome! If you encounter any issues or have suggestions for improvement, please feel free to open an issue or submit a pull request on the GitHub repository.

## License

This package is licensed under the [MIT License](https://opensource.org/licenses/MIT). See the LICENSE file for more details.
