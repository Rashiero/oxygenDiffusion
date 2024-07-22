import json
from PIL import Image,ImageSequence
import numpy as np
import sys

def main():
    # Scope
    fname = sys.argv[1]
    # Network
    with open(f"{fname}_final_network.json","r") as f:
        network = json.load(f)
    
    network['Diam'] = {}

    # Image
#    data = []
#    im = Image.open("Data/"+fname+".tif")
#    for i in ImageSequence.Iterator(im):
#        data.append(np.array(i))
#    del im
#    data = np.array(data)
#    data = data.astype("int8")

    # Define Routine to get radius
    for key in network['Vessel_list'].keys():
        network['Diam'][key] = 1e-4

    with open(f'{fname}_final_network.json','w') as f:
        json.dump(network,f)

if __name__ == "__main__":
    main()