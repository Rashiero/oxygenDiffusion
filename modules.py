import numpy as np
from parameters import *
import matplotlib.pyplot as plt
from scipy.fft import fftn, ifftn
import json

class Node():
    def __init__(self,*args):
        self.x = args[0] + margin_x//2
        self.y = args[1] + margin_y//2
        self.z = args[2] + margin_z//2
    def getCoord(self):
        return (self.x,self.y,self.z)
    def setBP(self,pressure:float):
        self.Pb = pressure
    def getBP(self):
        if "Pb" in self.__dict__.keys():
            return self.Pb
        else:
            raise ValueError("No value of Pressure assigned to the node")
    def setPressure(self,pressure:float):
        self.P = pressure
    def getPressure(self):
        if "P" in self.__dict__.keys():
            return self.P
        else:
            raise ValueError("No value of Pressure assigned to the node")
    def setFlow(self,flow):
        self.F = flow
    def getFlow(self):
        if "F" in self.__dict__.keys():
            return self.F
        else:
            print(self.getCoord())
            raise ValueError(f"No value of Flow assigned to the node")
    def __str__(self) -> str:
        return f"{self.getCoord()}"
    def __eq__(self,other):
        return self.getCoord() == other.getCoord()

class Vessel():
    def __init__(self,v1:Node = Node(0,0,0), v2:Node = Node(Nx-1,Ny-1,Nz-1)):
        # Start and Endpoints
        self.start = v1
        self.finish = v2
        self.fathers = []
        self.centerline = []
    def updateCenterline(self,center):
        for i in center:
            self.centerline.append(Node(*i))
    def __str__(self) -> str:
        return f"{self.start} -- {self.finish}"
    def getCoord(self):
        coords = []
        for i in self.centerline:
            coords.append(i.getCoord())
        return coords
    def setLength(self):
        length = 0
        vaso_old = self.start
        for node in self.centerline:
            d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
            ds = np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2)))
            length+=ds
            vaso_old = node
        node = self.finish
        d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
        ds = np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2)))
        length+=ds
        self.L = length    
    def setDiameter(self,diameter):
        self.D = diameter    
    def setResistance(self,res):
        self.R = res    
    def setFlux(self,Q):
        if Q < 0:
            self.start,self.finish = self.finish,self.start
            self.centerline = self.centerline[::-1]
            Q=-Q
        self.Q = Q   
    def getLength(self):
        return self.L    
    def getDiameter(self):
        return self.D    
    def getResistance(self):
        return self.R    
    def getFlux(self):
        return self.Q    
    def getStart(self):
        return self.start.getCoord()
    def getEnd(self):
        return self.finish.getCoord()

class Network():
    def __init__(self,path):
        with open(path,"r") as f:
            network = json.load(f)
        self.Vertex_list = network['Node_list']
        self.EdgePoints = network['Vessel_list']
        self.EdgeList = {}
        for vi,vf in self.EdgePoints.items():
            self.EdgeList[vi] = Vessel(Node(*self.Vertex_list[vf[0]]),Node(*self.Vertex_list[vf[1]]))
        for key in self.EdgeList.keys():
            self.EdgeList[key].updateCenterline(network['centerline'][key])
            self.EdgeList[key].setDiameter(network['Diam'][key])
    def setFathers(self):
        for eid in self.EdgeList.keys():
            for remaining_edges in self.EdgeList.keys():
                if self.EdgeList[eid].start == self.EdgeList[remaining_edges].finish:
                    self.EdgeList[eid].fathers.append(remaining_edges)
    def setLengths(self):
        for vaso in self.EdgeList.values():
            vaso.setLength()
    def setResistances(self):
        for vaso in self.EdgeList.values():
            vaso.setResistance(128*vaso.getLength()*mu/(np.pi*vaso.getDiameter()**4))
    def setFluxes(self,in_nodes = [], out_nodes = []):
        m_size = len(self.Vertex_list)
        map_key = {i:j for i,j in zip(self.Vertex_list.keys(),range(m_size))}
        adjacency_list = {x:[] for x in range(m_size)}
        for edge in self.EdgePoints.keys():
            start,end = self.EdgePoints[edge]
            adjacency_list[map_key[start]].append((end,edge))
            adjacency_list[map_key[end]].append((start,edge))
        b = np.zeros(m_size)
        M = np.zeros((m_size,m_size))
        for node in in_nodes:
            i = map_key[node]
            b[i] = PB0
            M[i][i] = 1
        for node in out_nodes:
            i = map_key[node]
            b[i] = PB1
            M[i][i] = 1
        for i in range(m_size):
            if b[i] == 0:
                for neighbors,edge in adjacency_list[i]:
                    M[i][map_key[neighbors]]+=1/self.EdgeList[edge].getResistance()
                    M[i][i]-=1/self.EdgeList[edge].getResistance()
        node_pressure = np.linalg.solve(M,b)
        for eid,nodes in self.EdgePoints.items():
            ini,fin = nodes
            self.EdgeList[eid].start.setBP(node_pressure[map_key[ini]])
            self.EdgeList[eid].finish.setBP(node_pressure[map_key[fin]])
        for key in self.EdgeList.keys():
            Q = (self.EdgeList[key].start.getBP() - self.EdgeList[key].finish.getBP())/self.EdgeList[key].getResistance()
            if Q < 0:
                vi,vf = self.EdgePoints[key]
                self.EdgePoints[key] = [vf,vi]
            self.EdgeList[key].setFlux(Q)


#####################
##### Functions #####
#####################

def updateConstants(new_val):
    global K
    global speed_const
    K = new_val
    speed_const = 1/K

def A_matrix():
    k1 = np.zeros((Nx+margin_x,Ny+margin_y,Nz+margin_z))
    k2 = np.zeros((Nx+margin_x,Ny+margin_y,Nz+margin_z))
    k3 = np.zeros((Nx+margin_x,Ny+margin_y,Nz+margin_z))
    k1[:,0,0] = np.fft.fftfreq(Nx+margin_x,d_x) * 2 * np.pi
    k2[0,:,0] = np.fft.fftfreq(Ny+margin_y,d_y) * 2 * np.pi
    k3[0,0,:] = np.fft.fftfreq(Nz+margin_z,d_z) * 2 * np.pi
    for i in range(Nz+margin_z):
        k2[0,:,i] = k2[0,:,0]
    for i in range(Nx+margin_x):
        k1[i,:,:] = k1[i,0,0]
        k2[i,:,:] = k2[0,:,:]
    k3[0,:,:] = k3[0,0,:]
    k3[:,:,:] = k3[0,:,:]
    k = 1/(D*(k1**2 + k2**2 + k3**2) + oxygen_degradation)
    return k

def flux(Pb,vaso:Vessel):
    return vaso.getFlux()*(Hd*Co*(Pb**n/(Pb**n+P50**n))+alfa*Pb)

def oxygenFlux(Po,net:Network):
    S = np.zeros((Nx+margin_x,Ny+margin_y,Nz+margin_z))
    for vessel in net.EdgeList.values():
        for point in vessel.centerline:
            indexes = point.getCoord()
            new_val = (point.getPressure() - Po[indexes])/K
            S[indexes] = new_val if new_val > 0 else 0
        new_val = (vessel.finish.getPressure() - Po[vessel.finish.getCoord()])/K
        S[vessel.finish.getCoord()] = new_val if new_val > 0 else 0
#            S[indexes] = speed_const*(point.getPressure() - Po[indexes])
#        S[vessel.finish.getCoord()] = speed_const*(vessel.finish.getPressure() - Po[vessel.finish.getCoord()])
    return S

def SteadyState(S,net:Network, A):
    s_hat = S.copy()
    for eid in net.EdgeList.keys():
        vaso = net.EdgeList[eid]
        vaso_old = vaso.start
        for node in vaso.centerline:
            d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
            s_hat[node.getCoord()] = S[node.getCoord()]*np.sqrt((d_step[0]*d_x)**2+(d_step[1]*d_y)**2+(d_step[2]*d_z)**2)/(d_x*d_y*d_z)
            vaso_old = node
        d_step = [y-x for x,y in zip(vaso_old.getCoord(),vaso.finish.getCoord())]
        s_hat[vaso.finish.getCoord()] = S[vaso.finish.getCoord()]*np.sqrt((d_step[0]*d_x)**2+(d_step[1]*d_y)**2+(d_step[2]*d_z)**2)/(d_x*d_y*d_z)
    s_hat = fftn(s_hat)
    pf_n = A*s_hat
    pf_n = np.real(ifftn(pf_n))
    return pf_n

def updateEntrance(net:Network):
    for key in net.EdgeList.keys():
        vaso = net.EdgeList[key]
#        print(f"{key} - {vaso.fathers} : Pressure: {vaso.start.getPressure()}")
        if np.isclose(vaso.getFlux(),0,atol=1e-20):
            vaso.start.setFlow(0.)
            continue
        n_fathers = len(vaso.fathers)
        if n_fathers == 0: # MUDAR
            vaso.start.setPressure(P0)
            vaso.start.setFlow(flux(P0,vaso))
        elif n_fathers == 1: # Tem caso?
            fid = vaso.fathers[0]
            vaso.start.setPressure(net.EdgeList[fid].finish.getPressure())
            vaso.start.setFlow(flux(vaso.start.getPressure(),vaso))
        else:
            vaso.start.setFlow(sum([net.EdgeList[i].finish.getFlow() for i in vaso.fathers]))
#            print(f"For {key}: Flow = {vaso.start.getFlow()} - Pressure should be between {flux(0,vaso)} and {flux(P0,vaso)}")
            vaso.start.setPressure(secantMethod(flux,0,P0,vaso,vaso.start.getFlow()))

def updateFlow(S,net:Network):#,Po):
    for vaso in net.EdgeList.values():
        vaso_old = vaso.start
        for node in vaso.centerline:
            d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
            ds = np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2)))
            node.F = vaso_old.getFlow() - S[vaso_old.getCoord()]*ds
#            if node.F < 0:
#                node.F = 0
#                node.setPressure(Po[node.getCoord()])
            node.F = node.F if node.F >=0 else 0
            vaso_old = node
        d_step = [y-x for x,y in zip(vaso_old.getCoord(),vaso.finish.getCoord())]
        ds = np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2)))
        vaso.finish.F = vaso_old.getFlow() - S[vaso_old.getCoord()]*ds
#        if vaso.finish.F < 0:
#            vaso.finish.F = 0
#            vaso.finish.setPressure(Po[vaso.finish.getCoord()])
        vaso.finish.F = vaso.finish.F if vaso.finish.F >=0 else 0

def setPb(net:Network,Po):
    for eid,Vaso in net.EdgeList.items():
        for node in Vaso.centerline:
            if node.getFlow() == 0:
                node.setPressure(Po[node.getCoord()])
            else:
                try:
                    node.setPressure(secantMethod(flux,0,P0,Vaso,val=node.getFlow()))
                except:
                    print(f"{eid} - {node}")
                    print(f"FAILED FOR FAILURE NODE WITH PRESSURE {node.getPressure()} AND FLOW {node.getFlow()}")
                    raise Exception("PRESSURE FAIL")
        if Vaso.finish.getFlow() == 0:
            Vaso.finish.setPressure(Po[Vaso.finish.getCoord()])
        else:
            Vaso.finish.setPressure(secantMethod(flux,0,P0,Vaso,val = Vaso.finish.getFlow()))

def secantMethod(func,x1,x2,params,val = 0, xacc = .6):
    MAX_ITER = 50
    fl = func(x1,params) - val
    f = func(x2,params) - val
    if(abs(fl) < abs(f)):
        rts = x1
        xl = x2
        fl,f=f,fl
    else:
        xl = x1
        rts = x2
    for _ in range(MAX_ITER):
        dx = (xl-rts)*f/(f-fl)
        xl = rts
        fl = f
        rts += dx
        f = func(rts,params) - val
        if abs(dx) < xacc or f == 0.0:
            return rts
    raise Exception(f"Maximum ({MAX_ITER}) number of iterations exceeded before Convergence")

def diagfig(Mat):
#    image = np.append(Mat[0,0,:],Mat[1,1,:])
#    for i in range(2,Mat.shape[0]):
#        image = np.append(image,Mat[i,i,:])
#    image = image.reshape(Mat.shape[0],Mat.shape[1])
#    return image
    return Mat[:,:,Nz//2]

def getFigure(Mat,flag = "show",title = "default"):
    fig = plt.figure()
    im = plt.imshow(diagfig(Mat))
    fig.colorbar(im,label = "Oxygen", orientation = "horizontal")
    if flag == "show":
        plt.show()
    elif flag == "save":
        plt.savefig(title)
    plt.close("all")

def generatePictures(Po,S,net,K,i):
    for eid,vasos in net.EdgeList.items():
        x_range = []
        vaso_old = vasos.start
        for node in vasos.centerline:
            d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
            x_range.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
            vaso_old = node
        d_step = [y-x for x,y in zip(vaso_old.getCoord(),vasos.finish.getCoord())]
        x_range.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
        x_range = np.array([0,*np.cumsum(x_range)])
        coord = [vasos.start.getCoord(),*vasos.getCoord(),vasos.finish.getCoord()]
        coord = (np.array([x[0] for x in coord]),np.array([x[1] for x in coord]),np.array([x[2] for x in coord]))
        # Po
        plt.subplot(2,2,1)
        plt.plot(x_range,Po[coord])
        plt.title(f"Po - Tissue Oxygen - {eid}")
        # Pb
        Pb = [vasos.start.getPressure()]
        Pb.extend([x.getPressure() for x in vasos.centerline])
        Pb.append(vasos.finish.getPressure())
        plt.subplot(2,2,2)
        plt.plot(x_range,Pb)
        plt.title(f"Pb - Blood Oxygen - {eid}")
        # S
        plt.subplot(2,2,3)
        plt.plot(x_range,S[coord])
        plt.title(f"S - Source - {eid}")
        # f
        flow = [vasos.start.getFlow(),*[x.getFlow() for x in vasos.centerline],vasos.finish.getFlow()]
        plt.subplot(2,2,4)
        plt.plot(x_range,flow)
        plt.title(f"f - Oxygen Flux - {eid}")
        #df/ds
        dfds = [0,*[(flow[i]-flow[i+1])/(x_range[i+1]-x_range[i]) for i in range(len(x_range)-1)]]
        plt.subplot(2,2,3)
        plt.plot(x_range,dfds)
        plt.savefig(f"Images/{eid}_{K}_{i}_plots.png")
    #    plt.show()
        plt.close("all")

def generatePictures_vessel(Po,S,net,K,i,eid):
    vasos = net.EdgeList[eid]
    x_range = []
    vaso_old = vasos.start
    for node in vasos.centerline:
        d_step = [y-x for x,y in zip(vaso_old.getCoord(),node.getCoord())]
        x_range.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
        vaso_old = node
    d_step = [y-x for x,y in zip(vaso_old.getCoord(),vasos.finish.getCoord())]
    x_range.append(np.sqrt(np.sum(np.power(np.array(d_step)*np.array([d_x,d_y,d_z]),2))))
    x_range = np.array([0,*np.cumsum(x_range)])
    coord = [vasos.start.getCoord(),*vasos.getCoord(),vasos.finish.getCoord()]
    coord = (np.array([x[0] for x in coord]),np.array([x[1] for x in coord]),np.array([x[2] for x in coord]))
    # Po
    plt.subplot(2,2,1)
    plt.plot(x_range,Po[coord])
    plt.title(f"Po - Tissue Oxygen - {eid}")
    # Pb
    Pb = [vasos.start.getPressure()]
    Pb.extend([x.getPressure() for x in vasos.centerline])
    Pb.append(vasos.finish.getPressure())
    plt.subplot(2,2,2)
    plt.plot(x_range,Pb)
    plt.title(f"Pb - Blood Oxygen - {eid}")
    # S
    plt.subplot(2,2,3)
    plt.plot(x_range,S[coord])
    plt.title(f"S - Source - {eid}")
    # f
    flow = [vasos.start.getFlow(),*[x.getFlow() for x in vasos.centerline],vasos.finish.getFlow()]
    plt.subplot(2,2,4)
    plt.plot(x_range,flow)
    plt.title(f"f - Oxygen Flux - {eid}")
    #df/ds
    dfds = [(flow[i]-flow[i+1])/(x_range[i+1]-x_range[i]) for i in range(len(x_range)-1)]
    dfds.append(dfds[-1])
    plt.subplot(2,2,3)
    plt.plot(x_range,dfds)
    plt.subplot(2,2,2)
    plt.plot(x_range,Po[coord])
    plt.subplot(2,2,4)
    plt.plot(x_range,[flux(P0,vasos)]*len(x_range))
    plt.savefig(f"Images/{eid}_{K}_{i}_plots.png")
#    plt.show()
    plt.close("all")

def save_vti_file(phi, nx, ny, nz, name):    
    pc_real = phi.real #+ psi.real
    with open(name + ".vti", "w" ) as my_file:
        my_file.write('<?xml version="1.0"?>')
        my_file.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        my_file.write('  <ImageData WholeExtent="0 '+str(nx)+' 0 '+str(ny)+' 0 '+str(nz)+'" Origin="0 0 0" Spacing ="1 1 1">\n')
        my_file.write('    <Piece Extent="0 '+str(nx)+' 0 '+str(ny)+' 0 '+str(nz)+'">\n') # dimensao da matriz x1 x2 y1 y2 z1 z2
        my_file.write('     <CellData>\n')
        my_file.write('     <DataArray Name="scalar_data" type="Float64" format="ascii">\n')
        my_file.write('     ')
        for iz in range(nz):
            pc_lista_novo = []
            for iy in range(ny):
                for ix in range(nx):
                    t = pc_real[ix,iy,iz]
                    pc_lista_novo.append(t)
            pc_string_novo = "    ".join([str(_) for _ in pc_lista_novo]) # criacao de uma string com os valores da lista
            my_file.write(pc_string_novo)
        my_file.write('\n         </DataArray>\n')
        my_file.write('      </CellData>\n')
        my_file.write('    </Piece>\n')
        my_file.write('</ImageData>\n')
        my_file.write('</VTKFile>\n')
        my_file.close()