import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def connected_components(adj_mat):
    nV = len(adj_mat)
    visited = [False for _ in range(nV)]
    cc = []
    for start_vertex in range(nV):
        if not visited[start_vertex]:
            stack = [start_vertex]
            comp = []
            while stack:
                v = stack.pop()
                if not visited[v]:
                    visited[v] = True
                    comp.append(v)
                    for i in range(nV): 
                        if adj_mat[v, i] == 1 and not visited[i]:
                            stack.append(i)
            cc.append(comp)
    return cc


with open('network_trial2.json','r') as f:
    network = json.load(f)

m_size = len(network['Node_list'])
key_list = np.array(list(network['Node_list'].keys()))
dist_m = np.zeros((m_size,m_size))

key_map = {key_list[i]:i for i in range(len(key_list))}
map_key = {j:i for i,j in key_map.items()}

perc_tot = m_size**2/2
perc_count = 0
for i in range(m_size):
    for j in range(i,m_size):
            dist_m[i,j] = np.sqrt(np.sum([(x-y)**2 for x,y in zip(network['Node_list'][key_list[i]],network['Node_list'][key_list[j]])]))
            dist_m[j,i] = dist_m[i,j]
            perc_count+=1
            print(f'{100*perc_count/perc_tot:.2f}%',end="\r")

teste = dist_m.copy()
teste[teste<2] = 1
for i in range(m_size):
    teste[i,i] = 0

teste[teste>=2] = 0

patterns = connected_components(teste)

for i in range(len(patterns)):
    if len(patterns[i]) > 1:
        edges = []
        for v in [map_key[x] for x in patterns[i]]:
            for key in network['Vessel_list'].keys():
                if v in network['Vessel_list'][key]:
                    edges.append(key)
        
        G = nx.from_pandas_edgelist(pd.DataFrame({'source':[network['Vessel_list'][x][0] for x in edges],
                                                  'target':[network['Vessel_list'][x][1] for x in edges]}))
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_labels(G, pos, font_size=8)
        nx.draw_networkx_edges(G, pos)
        labels = {tuple(sorted(network['Vessel_list'][i])):f"{i}={len(network['centerline'][i])}" for i in edges}
        nx.draw_networkx_edge_labels(G,pos,labels)
        plt.show()


adj_mat = np.zeros((m_size,m_size))
for key in network['Vessel_list'].keys():
    vi,vf = network['Vessel_list'][key]
    adj_mat[key_map[vf],key_map[vi]] = 1
    adj_mat[key_map[vi],key_map[vf]] = 1

count_degree = {i:0 for i in np.unique(np.sum(adj_mat,axis=0))}

for i in np.sum(adj_mat,axis=0):
    count_degree[i]+=1