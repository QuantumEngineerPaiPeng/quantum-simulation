# quantum-simulation
Codes to simulate spin dynamics in 1D spin chain

Except some early codes, most of the simulations are based on a matlab class "OperatorClass.m" in the "basic funcs" folder. The OperatorClass allows to store and operate (+-*/, diagonalize etc) the operators in symmetry sectors.
To project the operators onto symmetry sectors, the code will create and save the project operators in some .mat files and the path is written as the local path on my computer, which will lead to errors when run on a different computer. Simply modify the path variable will solve this error.
