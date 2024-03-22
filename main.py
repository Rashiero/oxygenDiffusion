from modules import *

##### Teste #####
error_dict = {}

##### Build World #####
VesselNetwork = Network("network.json")
Po = np.zeros((Nx,Ny,Nz))
S = np.zeros((Nx,Ny,Nz))
A = A_matrix()
# K = np.longfloat(1.2e10)
speed_const = 1/K

print("Network Built")

##### Blood Flow Properties #####
# Calculate Length of Centerlines
VesselNetwork.setLengths()

# Calculate Resistance
VesselNetwork.setResistances()

# Calculate Fluxes
VesselNetwork.setFluxes()
print("Flux calculated")

# Set initial conditions
for vessel in VesselNetwork.EdgeList.values():
    vessel.start.setPressure(P0)
    vessel.start.setFlow(flux(P0,vessel))
    for node in vessel.centerline:
        node.setPressure(P0)
        node.setFlow(flux(P0,vessel))
    vessel.finish.setPressure(P0)
    vessel.finish.setFlow(flux(P0,vessel))

VesselNetwork.setFathers()

# for i in VesselNetwork.EdgeList.keys():
#     print(f'{i} --- {VesselNetwork.EdgeList[i].fathers}')

##### Tissue Perfusion #####

j = 0
while K >= 3e8:
    print(f"K = {K:e} gamma = {speed_const:e}")
    err = 2e-3
    prev_err = 3e-3
    i = 0
    S = oxygenFlux(Po,VesselNetwork)
    error_dict[str(K)] = {"err":[],"alfafa":[]}

    while err > 1e-4 and not np.isclose(err/prev_err,1):
        alfafa = 0.5*100/(100+i)
        nS = S.copy()
        # Calculate Sources
        S = oxygenFlux(Po,VesselNetwork)
        S = nS + (S - nS)*alfafa
        nPo = Po.copy()
        # Calculate Steady State
        Po = SteadyState(S,VesselNetwork,A)
        # Update initial pressures
        updateEntrance(VesselNetwork)
        # Update Blood Oxygen Flow
        updateFlow(S,VesselNetwork)
        # Update Blood Oxygen Pressure
        setPb(VesselNetwork)
        
        # Prints and flow controls
        prev_err = err
        err = np.abs(nPo-Po).max()
        error_dict[str(K)]["err"].append(err)
        error_dict[str(K)]["alfafa"].append(alfafa)
        print(f"Iter == {i}, Erro == {err:.4e}, Atualização S == {np.abs(nS-S).max():.5e}", end="\r")
        if i%n_print == 0:
            getFigure(Po,"save",f"Images/{j}_{i}_O2field.png")
            generatePictures1(Po,S,VesselNetwork,j,i,'e21')
        i+=1
    print("")
    j+=1
    K-=2*10**(int(np.log10(K))-1)
    updateConstants(K)

save_vti_file(Po,Nx,Ny,Nz,"Tissue_Oxygen_Pressure")


with open("Images/error_dict.json",'w') as f:
    f.write(json.dumps(error_dict))