# Term Paper for Dynamic Programming

# Installing Kernel

> NOTE: It's only necessary to install dependencies if one doesn't have: pandas, scipy, numpy, numba etc.

Use pipenv to install dependencies:

- Navigate to the root directory of this folder. The **Pipfile** should be located at your path.
- Run the command: `pipenv install Pipfile` to install the python dependencies

Installing the kernel:

- Run following command in terminal: `pipenv run python -m ipykernel install --user --name=datascience --display-name "Python 3 (dypa-termpaper)`

# Running Code

We have attached an example notebook _example notebook_ that solves the model, and next simulates from the model. The model is solved for education level = 'college', and 1000 agents are simulated.

# Infrastructure of code

The 4 main files is placed in the folder _src_. These are:

- agent.py
- model.py
- parameters.py
- simulator.py

Where **parameters** contains the relevant parameter values, and **agent**, **model** and **simulator** contains classes that corresponds to their name.

Within the _src_ another folder is located called _modules_. This module, contains all relevant helper functions used in **agent**, **model**, **simulator** and **parameters**.

# Motivation

A study on the influence of financial knowledge on retirement wealth Inequality

Model inspired by:

- **Optimal Financial Knowledge and Wealth Inequality** by the authors _Annamaria Lusardi, Pierre-Carl Michaud, Olivia S. Mitchell_
- **Buffer Stock Model**
