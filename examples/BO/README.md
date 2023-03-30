# Bayesian Optimization Procedure

The Bayesian optimization (BO) procedure in this work was implemented using the BoTorch module built on PyTorch, which is a modern programming framework that enables flexible specification and optimization of probabilistic models [5]. 

![BO procedure](./figures/figure_S2.png)

## Initialization

For initialization, the design space for optimization should be defined, constructed by the bounds of the parameters, optimization target, and constraints. In this work, the variables are conditional and kinetic parameters involved in the elementary reactions in free radical copolymerization, where the value ranges were set based on the ranges reported in reference [6], as shown in the Table S2. 

After setting the design space, the initial dataset was established based on the Monte Carlo evaluation results of 28 randomly selected combinatorial parameters. Then, the Gaussian process model was trained by fitting the initial dataset to start the optimization loop. 

## Optimization Loop

After training the initial model, five sets of combinatorial parameters (candidates) to be evaluated by Monte Carlo can be outputted by optimizing the acquisition function for each iteration, where Expected Improvement [7] was employed as the acquisition function. Constraints were incorporated in the optimization process to exert bias on the selection of condition candidates. The Monte Carlo simulation results of the candidates (including total conversion, molecular weight (Mn), and molecular weight distribution (Ɖ)) were fed back to the original dataset to train an updated Gaussian process model for the next iteration. For a complete Bayesian optimization process, 200 iterations were conducted, giving a total of 1028 Monte Carlo evaluations (28 during initialization and 1000 during the optimization loop).
