# Datasets of the Gaussian model predictions of different FRP scenarios
This repository contains the Guassian model predictions results of four different polymerization methods optimized by Bayesian optimization. The filenames with "photo" or "therm" correspond to photo-initiated or thermally-initiated polymerization, while "disp" or "coup" correspond to disproportionation termination and coupling termination. The Gaussian models were trained based on the MC simulation results as described in **./datasets/BO_results**.

Each CSV file contains all the prediction results of the combinatorial parameter space for different FRP scenarios. Each row in the file represents a simulation, where the first 13 columns represent the input reaction parameters, and the last 3 columns represent the simulated molecular weight distribution, molecular weight, and total monomer conversion, respectively. The 13 reaction parameters are listed in the following order: log(*k*<sub>i</sub>), log(*k*<sub>paa</sub>), log(*k*<sub>pbb</sub>), log(*k*<sub>1</sub>), log(*k*<sub>2</sub>), log(*k*<sub>taa</sub>/*k*<sub>paa</sub>), log(*k*<sub>tbb</sub>/*k*<sub>pbb</sub>), log(*k*<sub>traa</sub>/*k*<sub>paa</sub>), log(*k*<sub>trab</sub>/*k*<sub>pab</sub>), log(*k*<sub>trba</sub>/*k*<sub>pba</sub>), log(*k*<sub>trbb</sub>/*k*<sub>pbb</sub>), log([M]/[Ini]), log([M]).
