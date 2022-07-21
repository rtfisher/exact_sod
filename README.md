#exact_sod by Robert Fisher

This program gives the exact analytical solution of the Sod shock tube test problem.

Users can use this program to test their hydrodynamics solvers against the exact analytical solution to check the accuracy of their solver.

It also contains an interactive plot library plotly https://plotly.com/python that helps to easily visulaize the data obtained. To run the program:

1. make sure you have a fortran compiler installed in your machine.
 run the following commands:
2.    ~ gfortran sod.f
3.    ~ ./a.out
4.    ~ python plot.py   

before executing the plot.py file, make sure you have plotly installed as well ( pip install plotly)

Have Fun!
