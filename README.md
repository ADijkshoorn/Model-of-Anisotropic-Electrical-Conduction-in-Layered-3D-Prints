# Model-of-Anisotropic-Electrical-Conduction-in-Layered-3D-Prints
Authors: A.P. Dijkshoorn, M. Schouten, S. Stramigioli, G.J.M. Krijnen 
Robotics and Mechatronics Group, Faculty of Electrical Engineering, Mathematics and Computer Science
University of Twente

Corresponding author: A.P. Dijkshoorn

Contact Information:
a.p.dijkshoorn@utwente.nl

University of Twente - Faculty of Electrical Engineering, Mathematics and Computer Science
P.O. Box 217
7500 AE Enschede
The Netherlands

***Introduction***
The code implements the model for anisotropic electrical conduction in layered 3D-prints, where it works for an arbitrary number of traxels (short for track-elements). It has the possibility to define the boundary conditions by means of a vector for the left boundaries and a vector for the right boundaries, where boundaries can be open, connected to other traxels, have input voltages and can be grounded. 
The code consists of various functions, where the "InitParameters" function allows for altering the material properties, parameters and boundary conditions.

From the parameters the code then calculates the voltages, current density and power dissipation for every traxel as well as the total resistance of the 3D-print. The model has been verified through simplified analytical examples and Finite Element Method simulations, showing good correspondence.
