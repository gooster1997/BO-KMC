# Demo test

This is a demo_test designed for quickly running one round of Bayesian optimization cycle. After installing the Monte Carlo program in the "**./Monte Carlo**" folder, users can run the "demo_test.py" file to explore the conditional space of FRP reactions under light-induced coupling termination. The script runs 20 cycles and the final output will be saved in the "output.csv" file, containing a total of 128 sets of simulation results (28 initialization sets + 20*5 optimization sets). The corresponding files have already been uploaded to this folder.

## Installation
To use this demo_test, you need to first compile the Monte Carlo program located in the "**./Monte Carlo**" folder. Please follow the instructions in the README file in the main folder to install and compile the program.

## Usage
After compiling the Monte Carlo program, you can run the "demo_test.py" file to start exploring the conditional space of FRP reactions under light-induced coupling termination. The script runs 20 cycles of Bayesian optimization and outputs the final results in the "output.csv" file.

To run the script, simply navigate to the directory where "demo_test.py" is located and type the following command in the terminal:
```
python demo_test.py
```
The script will start running and you will see the progress on the terminal. Once it's done, you can find the results in the "output.csv" file.
