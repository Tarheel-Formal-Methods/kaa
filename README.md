# Kaa
Kaa is a tool for reachability analysis of polynomial dynamical systems using parallelotope bundles.
It is a rewrite of the tool Sapo introduced by Dreossi, Dang, Piazza [(paper)](https://dl.acm.org/doi/abs/10.1145/2883817.2883838)
and formally through Dreossi's paper: 
[Sapo: Reachability Computation and Parameter Synthesis of Polynomial Dynamical Systems](https://arxiv.org/pdf/1607.02200.pdf)
Other papers elaborating on the techniques used in Kaa are:
[Parameter Synthesis using Parallelotopic Enclosure and Applications to Epidemic Models](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.707.7012&rep=rep1&type=pdf)

# Dependencies
Kaa relies only on the following python3 packages:

- [sympy](https://pypi.org/project/sympy/)
- [scipy](https://pypi.org/project/scipy/)
- [numpy](https://pypi.org/project/numpy/)
- [swiglpk](https://pypi.org/project/swiglpk/)

All of which can be installed through pip or through the package's corresponding page on PyPI.

# Running Sample Models and Examples.

The Juypter notebook named kaa-intro provides basic examples of the usage of kaa. The notebook contains all of the relevant code necessary to begin plotting the reachable sets and phase plots.

# Contents:
- [Summary of Files](md/explan.md)
