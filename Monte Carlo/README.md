# Monte Carlo Simulation Procedure

This README file describes the procedure of Monte Carlo simulation used in this work, which is based on the algorithm according to literature.1 Figure 1 describes the flow diagram of the Monte Carlo simulation in this work.
![flow diagram](./flow diagram.png)

## Initialization

The following steps were taken for initialization:
- Set the total volume V, molecular weights of monomers (MM1, MM2), rate constants of the related elementary reactions (kv, v = 1, 2, …, n)
- Set or calculate the numbers of the reaction reagents/products (initiator (nIni), monomers (nM1, nM2), radicals (nR1•, nR2•), dead polymer (nP)) based on the initial concentrations of initiators and monomers ([Ini]0, [M1]0 and [M2]0)
- Calculate the Monte Carlo rate constants (or microscopic rate constants) kv, MC based on their corresponding macroscopic rate constants kv by equation in the initialization section in Figure S1.

## Monte Carlo Simulation

The following steps were taken for Monte Carlo simulation:
- Calculate the reaction rates (Rv) and probabilities (Pv) to determine the time step τ and reaction channel μ for the next occurred reaction via the generation of two random numbers r1, r2
- For reactions involving the propagating radicals, one or two additional random numbers (r3 and/or r4) is generated to determine the specific chain (with distinct degree of polymerization and monomer sequence) for reaction
- Execute the reaction and update the numbers of the reactive species for the next iteration, until half consumption of the initiator was reached
- Calculate Mn, Ð and monomer conversions as the output.

## Elementary Reactions

The related elementary reactions used in this work are shown in Table S1.

## Reference
1. Gillespie, D. T., A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. J. Comput. Phys. 22, 403-434 (1976). 

Note: This README file is based on the description provided in the original work. For more details, please refer to the reference literature.
