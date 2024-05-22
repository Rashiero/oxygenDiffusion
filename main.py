from modules import *

##### Teste #####
error_dict = {}

##### Build World #####
VesselNetwork = Network("network.json")
Po = np.zeros((Nx,Ny,Nz))
S = np.zeros((Nx,Ny,Nz))
A = A_matrix()
# K = np.longfloat(1.2e10)
# speed_const = 1/K

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
try:
    while K >= 3e8:
        print(f"K = {K:e} gamma = {speed_const:e}")
        err = 2e-3
        prev_err = 3e-3
        i = 0
        S = oxygenFlux(Po,VesselNetwork)
        error_dict[str(K)] = {"err":[],"alfafa":[],"e21":{}}#,"e1":[],"e2":[],"e3":[],"e4":[],"e5":[],"e6":[],"e7":[],"e8":[],"e9":[],"e11":[],"e12":[],"e13":[],"e14":[],"e15":[],"e16":[],"e18":[],"e19":[],"e20":[],"e21":[],"e22":[],"e23":[],"e24":[],"e25":[]}
        while err > 1e-4:
            alfafa = 0.25#*100/(100+i)
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
            if np.isclose(err-prev_err,0):
                print("NOT CONVERGENT")
                break
            error_dict[str(K)]["err"].append(err)
            error_dict[str(K)]["alfafa"].append(alfafa)
            if i == 1:
                max_flow = flux(P0,VesselNetwork.EdgeList['e21'])
                ds = [0]
                vaso_old = VesselNetwork.EdgeList['e21'].start
                for node in VesselNetwork.EdgeList['e21'].centerline:
                    d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
                    ds.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
                    vaso_old = node
                d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
                ds.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
                coord = [VesselNetwork.EdgeList['e21'].start.getCoord(),*VesselNetwork.EdgeList['e21'].getCoord(),VesselNetwork.EdgeList['e21'].finish.getCoord()]
                coord = (np.array([x[0] for x in coord]),np.array([x[1] for x in coord]),np.array([x[2] for x in coord]))
                error_dict[str(K)]['e21']['Pb'] = [VesselNetwork.EdgeList['e21'].start.getPressure(),*[te.getPressure() for te in VesselNetwork.EdgeList['e21'].centerline],VesselNetwork.EdgeList['e21'].finish.getPressure()]
                error_dict[str(K)]['e21']['Po'] = [float(x) for x in Po[coord]]
                error_dict[str(K)]['e21']['S'] = [float(x) for x in S[coord]]
                error_dict[str(K)]['e21']['f'] = [VesselNetwork.EdgeList['e21'].start.getFlow(),*[te.getFlow() for te in VesselNetwork.EdgeList['e21'].centerline],VesselNetwork.EdgeList['e21'].finish.getFlow()]
                error_dict[str(K)]['e21']['ds'] = [float(x) for x in ds]
#                for key in VesselNetwork.EdgeList.keys():
#                    error_dict[str(K)][key] = [float((x-y)/K) for x,y in zip(,)]
#                    error_dict[str(K)][key] = [float((curr_flow-(x-y)*step/(K))/max_flow) for x,y,step,curr_flow in zip([VesselNetwork.EdgeList[key].start.getPressure(),*[te.getPressure() for te in VesselNetwork.EdgeList[key].centerline],VesselNetwork.EdgeList[key].finish.getPressure()],Po[coord],ds,[VesselNetwork.EdgeList[key].start.getFlow(),*[teee.getFlow() for teee in VesselNetwork.EdgeList[key].centerline],VesselNetwork.EdgeList[key].finish.getFlow()])]
            print(f"Iter == {i}, Erro == {err:.4e}, Atualização S == {np.abs(nS-S).max():.5e}\r", end = "")
            generatePictures_vessel(Po,S,VesselNetwork,j,i,'e21')
#        if i%n_print == 0:
#            getFigure(Po,"save",f"Images/{j}_{i}_O2field.png")
            i+=1
        print("")
        j+=1
#        K-=2*10**(int(np.log10(K))-1)
        K-=0.2*K
#    K-=0.2*K
        updateConstants(K)
except:
    pass
save_vti_file(Po,Nx,Ny,Nz,"Tissue_Oxygen_Pressure")


with open("Images/error_dict.json",'w') as f:
    f.write(json.dumps(error_dict))