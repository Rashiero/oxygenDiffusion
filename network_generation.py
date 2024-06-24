from PIL import Image,ImageSequence
import numpy as np
import json

flag = True

def neighborhood(data,ix,iy,iz,Nx,Ny,Nz,step = 1):
    neigh = []
    for ixi in range(max(0,ix-step),min(Nx,ix+step+1)):
        for iyi in range(max(0,iy-step),min(Ny,iy+step+1)):
            for izi in range(max(0,iz-step),min(Nz,iz+step+1)):
                if data[ixi,iyi,izi] != 0 and (ix != ixi or iy != iyi or iz != izi):
                    neigh.append((ixi,iyi,izi))
    return neigh

# data
data = []
fname = 'P8_Segmented'
im = Image.open(fname+"_skeleton.tif")
for i in ImageSequence.Iterator(im):
    data.append(np.array(i))

del im

data = np.array(data)
data = data.astype("int8")
Nx,Ny,Nz = data.shape


KNN = np.zeros((Nx,Ny,Nz))
perc_tot = Nx*Ny*Nz
perc_count = 0
for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            if data[ix,iy,iz] != 0:
                for ixi in range(max(0,ix-1),min(Nx,ix+2)):
                    for iyi in range(max(0,iy-1),min(Ny,iy+2)):
                        for izi in range(max(0,iz-1),min(Nz,iz+2)):
                            if ix != ixi or iy != iyi or iz != izi:
                                KNN[ix,iy,iz]+=data[ixi,iyi,izi]
            perc_count+=1
            print(f'{100*perc_count/perc_tot:.2f}%',end="\r")


network = {"Node_list":{},"Vessel_list":{},"centerline":{}}

current_vertex = 0
perc_tot = Nx*Ny*Nz
perc_count = 0
for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            if KNN[ix,iy,iz] != 0 and KNN[ix,iy,iz] != 2:
                current_vertex+=1
                key = 'v'+str(current_vertex)
                network["Node_list"][key]=[ix,iy,iz]
            perc_count+=1
            print(f'{100*perc_count/perc_tot:.2f}%',end="\r")

Search_space = KNN.copy()

queue = []
for V,coord in network["Node_list"].items():
    ix,iy,iz = coord
    if KNN[ix,iy,iz] == 1:
        queue.append(V)

current_edge = 0
checked_nodes = 0
full_nodes = len(network["Node_list"])
starting_vertex = queue[0]
while queue:
    v = queue.pop(0)
    checked_nodes+=1
    ix,iy,iz = network['Node_list'][v]
    Search_space[ix,iy,iz]=0
    neighbors = neighborhood(Search_space,ix,iy,iz,Nx,Ny,Nz)
    print(f'Checking node {v}; {checked_nodes} out of {full_nodes}')
    # Make edge
    for i in range(len(neighbors)):
        starting_vertex = neighbors[i]
        if Search_space[starting_vertex] != 0:
            key = 'e'+str(current_edge)
            centerline = []
            edge_size = 0
            while Search_space[starting_vertex] == 2:
                edge_size+=1
                print(f'neighbor {i} - size = {edge_size}',end="\r")
                centerline.append(starting_vertex)
                Search_space[starting_vertex] = 0
                ixi,iyi,izi = starting_vertex
                neigh2 = neighborhood(Search_space,ixi,iyi,izi,Nx,Ny,Nz)
                if len(neigh2) == 0:
                    flag = False
                    break
                starting_vertex = neigh2[0]
            # Find endpoint
            print('OUT')
            if flag:
                current_edge+=1
                endpoint = [x for x,y in network["Node_list"].items() if y == list(starting_vertex)][0]
                print(f'{v}:{endpoint} - neighbor {i} - total size = {edge_size+1}')
                network["Vessel_list"][key] = [v,endpoint]
                network["centerline"][key] = centerline
                queue.append(endpoint)
            else:
                flag = True



with open(f"{fname}_primary_network.json",'w') as f:
    json.dump(network,f)