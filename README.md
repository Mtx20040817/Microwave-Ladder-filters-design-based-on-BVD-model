# Microwave-Ladder-filters-design-based-on-BVD-model
This project aims to extract the parameter of the equivalent circuit based on the Generalized Low pass Cheybshev filter function with prescribed return loss and transmission zeros
Currently, this code only supports the design of symmetric filters with a start of series resonators. I will follow up the updated version regular.

Define_LP_function.m is used for defining the Generalized Low pass Cheybshev filter function

Parameter_extraction.m is used to extract the value of low pass FIR and admittance inverters.

LP_BP_transformation.m is used to do the low pass to bandpass transformation, which based on the narrow band approximation.

main.m is the part that calling these functions

Reference

[1] R. J. Cameron, C. M. Kudsia, and R. R. Mansour, Microwave Filters for Communication Systems: Fundamentals, Design and Applications. Hoboken, NJ, USA: Wiley, 2018. 

[2] A. Gimenez, J. Verdu, and P. De Paco Sanchez, “General synthesis methodology for the design of acoustic wave ladder filters and Duplexers,” IEEE Access, vol. 6, pp. 47969–47979, 2018. doi:10.1109/access.2018.2865808 


