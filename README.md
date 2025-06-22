# Phase-field codes to simulate crystal growth
The repository contains different implementations of phase-field model by Kobayashi (https://doi.org/10.1016/0167-2789(93)90120-P) to simulate crystal growth.

# C++/C_explicit
This folder contains C++ implementation of the model equations of the above reference. Both the phase-field and temperature equations are solved explicitly using finite volume method. This is a 2D serial code with flexibility to extend to 3D simulations.

Follow the steps to run the simulation code:

1. make clean

2. make

3. ./src/simulation input.in

Running the code will create a folder output where the output files are stored. The output files can be visualized using open-source packages such Paraview, gnuplot or Matplotlib.

# Python
This folder contains Python implementation of the model equation of the above reference. The governing equations are solved explicitly using a finite difference scheme. This is a 2D serial code.

Follow the steps to run the simulation code:

1. Install Conda :
    	Install Conda (Miniforge or Anaconda) so that the conda command is available in your shell.

2. Create a new environment :
        conda create -n phasefield_env python = 3.10
 
3. Activate the environment :
       conda activate phasefield

4. Install dependencies :
       pip install -r requirements.txt
     Note: Ensure the correct path to requirements.txt

5. Run the code :
     python3 src/main.py
     Note: Again ensure correct path to the src folder

Running the code will create a folder output where the output files and plots are stored

# C_explicit_implicit
This folder contains C implementation of the model equations of the above reference. The phase-field is solved explicitly and temperature equations are solved implicitly. The code is under construction and currently might not work as expected. 

Follow the steps to run the simulation code:

1. make clean

2. make

3. ./src/simulation input.in

Running the code will create a folder output where the output files are stored. The output files can be visualized using open-source packages such Paraview, gnuplot or Matplotlib.

