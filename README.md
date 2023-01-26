# **3D** model for **H**elium triplet **A**bsorption **T**ransmission **S**pectrum (3D HATS)

This is my Master's thesis code I've written under the supervision of Dr. Antonija Oklopcic. 

This code contains a 3D model of an exoplanetary atmosphere which is described by a 1D isothermal Parker atmospheric profile (for more details, see [Oklopcic & Hirata (2018)](https://iopscience.iop.org/article/10.3847/2041-8213/aaada9)). This 1D Parker profile is projected on a 3D spherical atmosphere, where day-to-night side winds and atmospheric rotation effects can be added. The schematic figure of the radiative transfer calculation is shown in the picture below. 

![Radiative transfer schematic](https://github.com/krsna336/3D_model_He_trip/blob/main/schematic_rt.png?raw=true)


Additionally, a limb darkening profile for the host star can be included. The final output from the code is the transmission spectrum of the exoplanet in Helium absoprtion triplet region (i.e., near 1083 nm).


**Files info** -
1. zeroth_3Dmodel.py -  Python file containing the 3D computational radiative transfer/atmopsheric model. 
2. test_parker_model.txt - A text file containing radial distributions of several parameters of WASP-69b (test planet). this will be called inside zeroth_3Dmodel.py
3. limb_darkening.py - A python code containing all the stellar limb darkening functions. This file is also called inside the main code zeroth_3Dmodel.py
