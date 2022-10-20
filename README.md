# 3D model for Helium triplet absorption transmission spectrum

This is my Master's thesis code I've written under the supervision of Dr. Antonija Oklopcic. 

This code contains a 3D model of an exoplanetary atmosphere which is described by a 1D isothermal Parker atmospheric profile (for more details, see [Oklopcic & Hirata (2018)](https://iopscience.iop.org/article/10.3847/2041-8213/aaada9)). This 1D Parker profile is projected on a 3D spherical atmosphere, where day-to-night side winds and atmospheric rotation effects can be added. The schematic figure of the radiative transfer calculation is shown in the picture below. 

![Radiative transfer schematic](https://github.com/krsna336/3D_model_He_trip/blob/main/schematic_rt.png?raw=true)


Additionally, a limb darkening profile for the host star can be included. The final output from the code is the transmission spectrum of the exoplanet in Helium absoprtion triplet region (i.e., near 1083 nm).
