# Bayesian Optimization Procedure

The Bayesian optimization (BO) procedure in this work was implemented using the BoTorch module built on PyTorch, which is a modern programming framework that enables flexible specification and optimization of probabilistic models.<sup>1</sup>

<img src="./flow-diagram.png" width="60%" height="60%">
Figure 1. Flow diagram of the BO framework integrated with Monte Carlo simulation.

## Initialization

For initialization, the design space for optimization should be defined, constructed by the bounds of the parameters, optimization target, and constraints. In this work, the variables are conditional and kinetic parameters involved in the elementary reactions in free radical copolymerization, where the value ranges were set based on the ranges reported in reference,<sup>2</sup> as shown in the Table 1. 

<img src="./value-ranges.png" width="60%" height="60%">
Table 1. Ranges of the evaluated reaction parameters of free radical copolymerization in during BO.

After setting the design space, the initial dataset was established based on the Monte Carlo evaluation results of 28 randomly selected combinatorial parameters. Then, the Gaussian process model was trained by fitting the initial dataset to start the optimization loop. 

## Optimization Loop

After training the initial model, five sets of combinatorial parameters (candidates) to be evaluated by Monte Carlo can be outputted by optimizing the acquisition function for each iteration, where Expected Improvement was employed as the acquisition function. Constraints were incorporated in the optimization process to exert bias on the selection of condition candidates. The Monte Carlo simulation results of the candidates (including total conversion, molecular weight (*M*<sub>n</sub>), and molecular weight distribution (*Ɖ*)) were fed back to the original dataset to train an updated Gaussian process model for the next iteration. For a complete Bayesian optimization process, 200 iterations were conducted, giving a total of 1028 Monte Carlo evaluations (28 during initialization and 1000 during the optimization loop).

## References
1. Balandat, M., et al., BoTorch: a framework for efficient Monte-Carlo Bayesian optimization. Proc. Adv. Neural Inf. Process. Syst. 33, 21524-21538 (2020).
2. Brandrup, J., et al., Polymer handbook. (Wiley New York, 1999).
