# Modelling Tissue Irrigation in the Neonatal Mouse Brain Cortex

The following code was developed in the context of the master in computational biology for the simulation of oxygen diffusion in the brain. The code was developed by Rafael Berton Correia Ramos under the guidance of [Rui Travasso](https://biologicalmodelling.com/people/) and [Vanessa Coelho-Santos](https://www.vcoelhosantoslab.com/team). The entire process takes as input a binary segment image of brain vasculature, process the images, extract vessel information and then simulates diffusion.
The codes developed are ment to be used in sequence and to aid in the process there is the `pipeline.sh` Bash that sequentially calls them.

**TODO**
**Important**
This code is not generalized. It makes the assumption of the input data being the ones in the `Data/` folder. As such it already assumes the names and the sizes. As the images are Z-stacks with differing Z-axis, part of the `pipeline.sh` involves calling `sed` to alter the parameter files accordingly to the name of the image used. Generalizing to accept any input will revolve mainly in the creation of variables where the names and sizes of the images are called and a function that would extract the sizes from the original image and input it into the parameters file. However this is not the scope of this project, as its first aim was to simulate blood flow to this specified data. 

# Files

## Image processing steps
### image_processing.py 

This routine processes the vessel network, cleaning the image and extracting the center of vessels via a skeletonization algorithm

### network_generation.py

This routine utilizes the centerline image to extract the network of blood vessels. It is stored in a .json file used to generate the python object.

### network_cleanup.py 

This routine processes the network into a compatible topology, reshaping all connections into simple bifurcations such that no node has higher degree that 3. It also trims smaller vessels that are artifacts from the skeletonization algorithm.

### calculate_diameter.py

This routine is not yet implemented. As the diameter is necessary for the `Network` class construction, this routine currently just updates the vessel network with a homogeneous sensible diameter value.

## Simulation

### parameters.py

This file calls a few necessary libraries and define parameters used on simulation. The rationale behind separating this script is to make parameters search easier, without altering the source code, as well as concentrating the information necessary to replicate a run of the program.

### modules.py

 This contains the bulk of the logic. The vessel network information is handled by the `Network` class, which in turn calls the `Node` and `Vessel` classes.

The tissue is best modelled as a 3D field, and is represented as a matrix. As consequence, the functions not encapsulated by the classes tend to either be related to tissue or to the crosstalk between vasculature and tissue.

### main.py 

This is the algorithm that utilizes the network information to simulate oxygen supply from diffusion. It imports `modules.py`, a file containing the definition of classes and functions. These modules in turn import the `parameters.py` file, containing the physical constants and other parameters.

The simulation presents some instability regarding the vessel wall permeability to Oxygen. The main loop first stabilises the oxygen field for the first iteration and then proceeds to the main simulation, where the K parameter will be decreased incrementally and in a adaptive step fashion.

## Others

### calculate_diameter.py
**TODO**
The objective of this routine is to use the vessel location information and extract from the original image the average radius or diameter of the vessel. The basic strategy consists of finding the volume of the vessel segment and in possession of the length of the vessel, fit into a cylinder volume formula to extract the radius 

### see_entrypoints.py

This algorithm isn't called, but it can aid in identifying the entry points, as it generates a image of the skeleton network with all possible incoming/outgoing vessels entries/exits enlarged. This also has a subroutine that extracts the `Nx`, `Ny` and `Nz` of each image. Here the caveat applies that the images are called by name and are not generalized.