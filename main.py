from modules import *
import copy

##### Teste #####
error_dict = {}
alfafa = 0.25
fname = "P8_Segmented"
##### Build World #####
VesselNetwork = Network(f"{fname}_final_network.json")
Po = np.zeros((Nx,Ny,Nz))
S = np.zeros((Nx,Ny,Nz))
A = A_matrix()

print("Network Built")

##### Blood Flow Properties #####
# Calculate Length of Centerlines
VesselNetwork.setLengths()

# Calculate Resistance
VesselNetwork.setResistances()

# Calculate Fluxes

in_nodes = []
out_nodes = []
for key in VesselNetwork.Vertex_list.keys():
    if VesselNetwork.Vertex_list[key][0] < 2:
        in_nodes.append(key)
    elif VesselNetwork.Vertex_list[key][0] >= 130:
        out_nodes.append(key)

degrees = {x:0 for x in VesselNetwork.Vertex_list.keys()}
for key in VesselNetwork.EdgePoints.keys():
    vi,vf=VesselNetwork.EdgePoints[key]
    degrees[vi]+=1
    degrees[vf]+=1

one_degree = [x for x,i in degrees.items() if i==1]
in_nodes = [x for x in in_nodes if x in one_degree]
out_nodes = [x for x in out_nodes if x in one_degree]

VesselNetwork.setFluxes(in_nodes,out_nodes)
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
err = 2e-3
prev_err = 3e-3
i = 0
error_dict[str(K)] = {"err":[],"alfafa":[]}
print(f"K = {K:e} gamma = {speed_const:e}")
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
    setPb(VesselNetwork,Po)
        
    # Prints and flow controls
    prev_err = err
    err = np.abs(nPo-Po).max()
    if np.isclose(err-prev_err,0):
        print("NOT CONVERGENT")
        break
    error_dict[str(K)]["err"].append(err)
    error_dict[str(K)]["alfafa"].append(alfafa)
    print(f"Iter == {i}, Erro == {err:.4e}, Atualização S == {np.abs(nS-S).max():.5e}\r", end = "")
    i+=1

print("Stabilized for first iteration")

j = 0
while K >= 3e8:
    step_flag = True
    step = 0.2*K
#    step = 2*10**(int(np.log10(K))-1)
    while step_flag:
        K1 = K-step
        updateConstants(K1)
        VN1 = copy.deepcopy(VesselNetwork)
        Pot = Po.copy()
        St = oxygenFlux(Pot,VN1)
        err1 = 0
        err = 2e-3
        prev_err = 3e-3
        i=0
        while err > 1e-4 and err1<=1.1:
            nS = St.copy()
            St = oxygenFlux(Pot,VN1)
            St = nS + (St - nS)*alfafa
            nPo = Pot.copy()
            Pot = SteadyState(St,VN1,A)
            updateEntrance(VN1)
            updateFlow(St,VN1)
            setPb(VN1,Pot)
            prev_err = err
            err = np.abs(nPo-Pot).max()
            if np.isclose(err-prev_err,0):
                print("NOT CONVERGENT")
                break
            for vaso in VN1.EdgeList.values():
                max_flux = flux(P0,vaso)
                vaso_flux = [vaso.start.getFlow()/max_flux,*[x.getFlow()/max_flux for x in vaso.centerline],vaso.finish.getFlow()/max_flux]
                err1 = max(vaso_flux) if max(vaso_flux) > err1 else err1
            print(f"Iter == {i}, Erro == {err:.4e}, Atualização S == {np.abs(nS-S).max():.5e}\r", end = "")
            i+=1

#        print(f"K = {K:.4e}\tK1 = {K1:.4e}\tstep = {step:.4e}\terr = {err1:.4e}\r", end = "")
        if err1 > 1.1:
#            step-=(1-err1)/err1*step
            step/=err1
        else:
            step_flag = False
    print("")
    K-=step
    updateConstants(K)
    print(f"K = {K:e} gamma = {speed_const:e}")
    err = 2e-3
    prev_err = 3e-3
    i = 0
    S = oxygenFlux(Po,VesselNetwork)
    while err > 1e-4:
        nS = S.copy()
        S = oxygenFlux(Po,VesselNetwork)
        S = nS + (S - nS)*alfafa
        nPo = Po.copy()
        Po = SteadyState(S,VesselNetwork,A)
        updateEntrance(VesselNetwork)
        updateFlow(S,VesselNetwork)
        setPb(VesselNetwork,Po)
        prev_err = err
        err = np.abs(nPo-Po).max()
        if np.isclose(err-prev_err,0):
            print("NOT CONVERGENT")
            break
        print(f"Iter == {i}, Erro == {err:.4e}, Atualização S == {np.abs(nS-S).max():.5e}\r", end = "")
        i+=1
    j+=1
    print("")


save_vti_file(Po,Nx,Ny,Nz,"Tissue_Oxygen_Pressure")