Description:

This program is designed to completely automate the process of two-loop self energy calculation taking an input of a FeynArts model file all the way through to a numerical result for the self energy.

We make use of the following tools to complete this process:

FeynArts (hep-ph/0012260) is used to generate the two-loop amplitudes
FeynCalc (ArXiv:1601.01167) is used to reduce the amplitudes
TARCER (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
TSIL (hep-ph/0501132) is used to evaluate the basis integrals

Required programmes:

To use this tool you need to have:

Mathematica
FeynCalc including a FeynArts installation (ships with by default)  # add more details here
TSIL


Instructions:

See user manual in documentation/mass_builder.pdf for installation and usage instructions