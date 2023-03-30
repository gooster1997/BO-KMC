# Bayesian Optimization for Free Radical Polymerization (BO-FRP)

This is a framework for exploring the condition space of free radical polymerization (FRP) using Bayesian optimization and Monte Carlo simulation. It can be used to find optimal polymerization conditions that result in a narrow molecular weight distribution. The Bayesian optimization part of the framework is developed based on the botorch module, which can be found at https://github.com/pytorch/botorch. For more details, please see our article which is currently under submission (the link will be added later).

## Requirements
The following modules are required to run the program:
- python==3.7.0
- numpy==1.21.5
- pandas==1.3.4
- torch==1.10.2+cu113
- torchvision==0.11.3+cu113
- gpytorch==1.6.0
- botorch==0.6.0

## Usage
1. (Installation)Navigate to the **./Monte Carlo** directory which contains Monte Carlo simulation algorithms for different types of FRP reactions compiled with the C language. In the filename of these algorithm files, "photo" or "therm" represents photoinduced or thermal initiation, and "disp" or "coup" represents disproportionation termination or coupling termination. Run the command ```python setup.py install``` to import these algorithm files into the Python environment for calling. It will take only seconds for the installation.
2. In the **./examples/BO** directory, Bayesian optimization .py files are available for optimizing suitable polymerization conditions for a specified type of FRP reaction to achieve any synthesis goal, including narrow molecular weight distribution. You can also modify the code to achieve other optimization goals and constraints. Specify the constrained objective in the following statement:
```
constrained_obj = ConstrainedMCObjective(
    objective=obj_callable,
    constraints=[constraint1_callable, constraint2_callable],
)
```
Additionally, the condition space on which optimization is based can also be arbitrarily specified. You can modify it in the following statement:
```
bounds = torch.tensor([[-8,-2,-2,-2,-2,1,1,-5,-5,-5,-5,-3,-1], [-6,5,5,2,2,7,7,-2,-1,-1,-2,-1,1]], device=device, dtype=dtype)
```
After reaching the specified number of iterations, optimization will stop and the Monte Carlo simulation results obtained during optimization will be output as a .csv file in the current folder.

3. In the **./examples/prediction** directory, the prediction.py file can be used to establish a Gaussian model based on the Monte Carlo simulation results for predicting the entire combination condition space. Prepare the "dataset.csv" file in the directory with the same format as the output file of Bayesian optimization. The final program will output the predicted values for all conditions corresponding to the set combination condition space and save the results as a "prediction_results.csv" file. You are free to set the structure (range, interval, etc.) of the combination condition space.

## Dataset
The Bayesian optimization and Gaussian model prediction results obtained in this article are stored in the **./datasets** folder.

## Contact
If you have any questions about this program, please contact chenmao@fudan.edu.cn or ygu19@fudan.edu.cn.



