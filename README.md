# Modelling blood flow in brain capillaries

The following codes are utilized in the processing of microscopy images and simulation of oxygen distribution in tissue. They are ment to be used in sequence and to aid in the process there is the `pipeline.sh` Bash that sequentially calls them.

**TODO**
**Important**
This code is not generalized. It makes the assumption of the input data being the ones in the `Data/` folder. As such it already assumes the names and the sizes. As the images are Z-stacks with differing Z-axis, part of the `pipeline.sh` involves calling `sed` to alter the parameter files accordingly to the name of the image used. Generalizing to accept any input will revolve mainly in the creation of variables where the names and sizes of the images are called and a function that would extract the sizes from the original image and input it into the parameters file. However this is not the scope of this project, as its first aim was to simulate blood flow to this specified data. 

# Files

## image_processing.py 

This routine processes the vessel network, cleaning the image and extracting the center of vessels via a skeletonization algorithm

## network_generation.py

This routine utilizes the centerline image to extract the network of blood vessels. It is stored in a .json file used to generate the python object.

## network_cleanup.py 

This routine processes the network into a compatible topology, reshaping all connections into simple bifurcations such that no node has higher degree that 3. It also trims smaller vessels that are artifacts from the skeletonization algorithm.

## calculate_diameter.py

**TODO**
This routine is not yet implemented. Its objective is to use the vessel location information and extract from the original image the average radius or diameter of the vessel. As the diameter is necessary for the `Network` class construction, this routine currently just updates the vessel network with a homogeneous sensible diameter value.

## main.py 

This is the algorithm that utilizes the network information to simulate oxygen supply from diffusion. It imports `modules.py`, a file containing the definition of classes and functions. These modules in turn import the `parameters.py` file, containing the physical constants and other parameters.

## see_entrypoints.py

This algorithm isn't called, but it can aid in identifying the entrypoints, as it generates a image of the skeleton network with all possible incoming/outgoing vessels entries/exits enlarged. This also has a subroutine that extracts the `Nx`, `Ny` and `Nz` of each image. Here the caveat applies that the images are called by name and are not generalized.