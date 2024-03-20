from modules import *
import json

##### Teste #####
error_dict = {}

##### Build Vessel Network #####
with open("network.json","r") as f:
    teste = json.load(f)
VesselNetwork = Network(teste['Node_list'],teste['Vessel_list'],teste['centerline'])
Po = np.zeros((Nx,Ny,Nz))
S = np.zeros((Nx,Ny,Nz))

print("Network Built")
# Get Vessel Diameters
i = 0
Diam = teste['Diam']
for vaso in VesselNetwork.EdgeList.values():
    vaso.setDiameter(Diam[i])
    i+=1
##### Blood Flow Properties #####
# Calculate Length of Centerlines
for vaso in VesselNetwork.EdgeList.values():
    vaso.setLength()
# Calculate Resistance
for vaso in VesselNetwork.EdgeList.values():
    vaso.setResistance(128*vaso.getLength()*mu/(np.pi*vaso.getDiameter()**4))
# Create adjacency list
node_codex = {}
index = 0
for i in teste['Node_list'].keys():
    node_codex[i]=index
    index+=1

adjacency_list = {x:[] for x in range(len(teste['Node_list']))}
for edge in VesselNetwork.EdgePoints.keys():
    start,end = VesselNetwork.EdgePoints[edge]
    adjacency_list[node_codex[start]].append((end,edge))
    adjacency_list[node_codex[end]].append((start,edge))
# Build system
# Kirchoff Current law
b = [0]*len(teste['Node_list'])
b[0] = PB0
b[-1] = PB1
# Flow Matrix
M = [[0]*len(teste['Node_list']) for _ in range(len(teste['Node_list']))]
M[0][0] = 1
M[-1][-1] = 1
for i in range(1,len(teste['Node_list'])-1):
    for neighbors,edge in adjacency_list[i]:
        M[i][node_codex[neighbors]]+=1/VesselNetwork.EdgeList[edge].getResistance()
        M[i][i]-=1/VesselNetwork.EdgeList[edge].getResistance()
# Solve system
node_pressure = np.linalg.solve(M,b)
for eid,nodes in VesselNetwork.EdgePoints.items():
    ini,fin = nodes
    VesselNetwork.EdgeList[eid].start.setBP(node_pressure[node_codex[ini]])
    VesselNetwork.EdgeList[eid].finish.setBP(node_pressure[node_codex[fin]])
# Calculate Volumetric Blood Flow via Ohm's Law
for vaso in VesselNetwork.EdgeList.values():
    vaso.setFlux((vaso.start.getBP() - vaso.finish.getBP())/vaso.getResistance())
print("Flux calculated")

for vessel in VesselNetwork.EdgeList.values():
    vessel.start.setPressure(P0)
    vessel.start.setFlow(flux(vessel.start.getPressure(),vessel))
    for node in vessel.centerline:
        node.setPressure(P0)
    vessel.finish.setPressure(P0)

edge_order = ["e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12","e13","e14","e15","e16","e17","e18","e19","e21","e20","e23","e24","e25"]
VesselNetwork.setFathers()
for i in VesselNetwork.EdgeList.keys():
    print(f'{i} --- {VesselNetwork.EdgeList[i].fathers}')

##### Tissue Perfusion #####
# Create Grids
Po = np.zeros((Nx,Ny,Nz))
A = A_matrix()
K = np.longfloat(1.2e10)
speed_const = 1/K
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
        Po = SteadyState(S,VesselNetwork,edge_order,A)
        # Update Blood Oxygen Flow
        updateFlow(S,VesselNetwork,edge_order)
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
            generatePictures(Po,S,VesselNetwork,j,i)
        i+=1
    print("")
    j+=1
    K-=2*10**(int(np.log10(K))-1)
    speed_const = 1/K
save_vti_file(Po,Nx,Ny,Nz,"Tissue_Oxygen_Pressure")


with open("Images/error_dict.json",'w') as f:
    f.write(json.dumps(error_dict))    