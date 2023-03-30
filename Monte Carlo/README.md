# Monte Carlo Simulation Procedure

This README file describes the procedure of Monte Carlo simulation used in this work, which is based on the algorithm according to literature.<sup>1</sup> Figure 1 describes the flow diagram of the Monte Carlo simulation in this work.

<img src="./flow-diagram.png" width="60%" height="60%">
Figure 1. Flow diagram of the algorithm designed for the kinetic Monte Carlo simulations of free radical copolymerizations.

## Initialization

The following steps were taken for initialization:
- Set the total volume V, molecular weights of monomers (M<sub>M1</sub>, M<sub>M2</sub>), rate constants of the related elementary reactions (*k*<sub>v</sub>, v = 1, 2, …, n)
- Set or calculate the numbers of the reaction reagents/products (initiator (n<sub>Ini</sub>), monomers (n<sub>M1</sub>, n<sub>M2</sub>), radicals (n<sub>R1•</sub>, n<sub>R2•</sub>), dead polymer (n<sub>P</sub>)) based on the initial concentrations of initiators and monomers ([Ini]<sub>0</sub>, [M1]<sub>0</sub> and [M2]<sub>0</sub>)
- Calculate the Monte Carlo rate constants (or microscopic rate constants) *k*<sub>v</sub>, MC based on their corresponding macroscopic rate constants *k*<sub>v</sub> by equation in the initialization section in Figure 1.

## Monte Carlo Simulation

The following steps were taken for Monte Carlo simulation:
- Calculate the reaction rates (*R*<sub>v</sub>) and probabilities (*P*<sub>v</sub>) to determine the time step τ and reaction channel μ for the next occurred reaction via the generation of two random numbers *r*<sub>1</sub>, *r*<sub>2</sub>
- For reactions involving the propagating radicals, one or two additional random numbers (*r*<sub>3</sub> and/or *r*<sub>4</sub>) is generated to determine the specific chain (with distinct degree of polymerization and monomer sequence) for reaction
- Execute the reaction and update the numbers of the reactive species for the next iteration, until half consumption of the initiator was reached
- Calculate *M*<sub>n</sub>, *Ð* and monomer conversions as the output.

## Elementary Reactions

The related elementary reactions used in this work are shown in Table 1.

<img src="./elementary-reactions.png" width="60%" height="60%">
Table 1. Elementary reactions considered in the kinetic Monte Carlo simulation of thermal-induced free radical copolymerization.

## Reference
1. Gillespie, D. T., A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. J. Comput. Phys. 22, 403-434 (1976). 

Note: This README file is based on the description provided in the original work. For more details, please refer to the reference literature.
