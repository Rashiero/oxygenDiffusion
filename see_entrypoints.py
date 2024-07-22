from PIL import Image,ImageSequence
import numpy as np
import json
import skimage
import sys

fname = sys.argv[1]
Nx,Ny,Nz = 133,512,512

with open(f"{fname}_final_network.json","r") as f:
    network = json.load(f)

# Get degrees
degs = {v:0 for v in network['Node_list'].keys()}
for key in network['Vessel_list'].keys():
    vi,vf = network['Vessel_list'][key]
    degs[vi]+=1
    degs[vf]+=1

# Set Degrees 1 as point
im = np.zeros((Nx,Ny,Nz))
for key in [x for x in degs.keys() if degs[x] == 1]:
    im[tuple(network['Node_list'][key])] = 1

# Diffusion
im = skimage.morphology.binary_dilation(im)
im = skimage.morphology.binary_dilation(im)
im = skimage.morphology.binary_dilation(im)

# Add centerlines
for key in network['Vessel_list'].keys():
    vi,vf = network['Vessel_list'][key]
    im[tuple(network['Node_list'][vi])] = 1
    im[tuple(network['Node_list'][vf])] = 1
    for i in network['centerline'][key]:
        im[tuple(i)] = 1

skimage.io.imsave(f"{fname}_entrypoints.tiff",im.astype("uint8"))

def dims():
    with open("Dimensions.txt","w") as f:
        f.write(f"File\tx\ty\tz\n")
        for fname in ["P8_Segmented","P9_Segmented","P10_Segmented","P11_Segmented","P12_Segmented"]:
            data=[]
            im = Image.open("Data/"+fname+".tif")
            for i in ImageSequence.Iterator(im):
                data.append(np.array(i))
            del im
            data = np.array(data)
            data = data.astype("int8")
            Nx,Ny,Nz = data.shape
            f.write(f"{fname}\t{Nx}\t{Ny}\t{Nz}\n")