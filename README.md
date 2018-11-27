# measuring.coevolution
Sample scripts and data to demonstrate the measurement of coevolution. A tutorial for this method is available [here](https://bobweek.github.io/measuring_coevolution.html).

  - __toju_pairs.Rda:__ data set containing pairs of local mean traits from multiple populations.
  - __toju_pars.Rda:__ data set containing background parameters of model.
  
Repository also includes original scripts used to produce figures and results of the associated manuscript:

  - __estimate.R:__ script used to estimate strengths of biotic selection and optimal offsets for each system analyzed.
  - __functions.R:__ script contains equilibrium expressions and utility functions used in other scripts.
  - __power_sample.R:__ script generates statistical power of our method as a function of sample size.
  - __power_strength.R:__ script generates statistical power of our method as a function of strength of coevolution.
  - __type1_sample.R:__ script generates type-1 error rate of our method as a function of sample size.
  - __type1_sample.R:__ script generates type-1 error rate of our method as a function of unilateral selection strength.
  - __regression_sample.R:__ script generates linear regression statistics of our method as a function of sample size.
  - __nn__: folder contains scripts used to analyze perfomance of method under non-normal data.
  - __gf__: folder contains scripts used to analyze perfomance of method in the presence of gene-flow.
  - __err__: folder contains scripts used to analyze perfomance of method with error in estimates of abiotic optima.
