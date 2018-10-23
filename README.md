# GAMMA_RAY_CMFGEN

This code is integrated into CMFGEN and is to create a synthetic gamma-ray spectrum and energy energy deposition in supernovae (thermonuclear and core-collapse). Spherical symmetry is assumed. 

## Running the Code

The code is wrirten in fortran 90/95 formatting. Makefiles need to be modified to reflect the users directory (Makefile\_definitions). Original code is compiled using the PGF 95 compiler. Gfortran will work with some massaging. 

Once the code is compiled, you can run the test model from the command line. Go into the test model directory and run the code using "batch.sh &". However, you may want to change certain parameters the gamma-ray part specifically uses to run (see file GAMRAY\_PARAMS). 
