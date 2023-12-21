# Bayesian Optimization for kinetic Monte Carlo simulation of Free Radical Polymerization (BO-KMC)

This project is a space exploration framework that integrates Bayesian optimization and dynamic Monte Carlo simulation. It aims to find suitable sets of reaction parameters for free radical polymerization to achieve lower dispersity.

## Installation Steps
## 1.Install Dependencies
This project is specifically designed to run in a Linux environment. Installing the KMC program on a Windows environment might result in errors.
Ensure that the following dependencies and their versions are installed in your environment:
- python==3.11
- numpy==1.26.0
- pandas==2.0.3
- pytorch=2.0.1
- gpytorch==1.6.0
- botorch==0.6.0 (for detailed installation guide on botorch, please see https://github.com/pytorch/botorch)

## 2.Clone the Project
```
git clone https://github.com/gooster1997/BO-KMC.git
cd BO-KMC
```
## 3.Install KMC Functions
Before running the main program, it's necessary to install the KMC calculation files using pypa/build. Follow these steps:
```
python -m build
cd dist
pip install YOUR_FILE_NAME.whl
```
Replace YOUR_FILE_NAME.whl with the actual name of the generated .whl file in the dist directory.
Afterward, the installation should be completed.

## 4.Running the Project
Run the BO_demo.ipynb file to execute this project.

## Contribution
This project credits Yu Gu, Tianyi Gao, and Mao Chen from the Department of Macromolecules, Fudan University. If you have any suggestions or find a bug, feel free to raise an issue or contact ygu19@fudan.edu.cn or chenmao@fudan.edu.cn.



