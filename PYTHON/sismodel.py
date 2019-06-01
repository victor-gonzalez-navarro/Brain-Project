import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import os


# ------------------------------------------------------------------------------------------------- READ NETWORK
graph = 'Coactivation.net'
print('Computing the network: '+graph)
G = nx.Graph(nx.read_pajek('NetworksPajek/' + graph))

# ------------------------------------------------------------------------------------------------- DRAW NETWORK
# fig = plt.figure()
# nx.draw_networkx(G, with_labels=False, node_size=20)
# limits = plt.axis('off')
# plt.suptitle('Network: ' + graph + '\n\n ', fontsize=16)
# fig.savefig('./ImagesSiS/' + graph + '.png')

# ------------------------------------------------------------------------------------------------ SIS ALGORITHM
Nnodes = len(G.nodes)
bet = 0.06
posible_betas = [bet]
ro_init = 0.2
Tmax = 500
Ttrans = Tmax-100

maxw = 0
for key1 in G._adj:
    for key2 in G._adj[key1]:
        if G._adj[key1][key2]['weight'] > maxw:
            maxw = G._adj[key1][key2]['weight']
posible_betas[0] = posible_betas[0] / maxw


fig = plt.figure(figsize=(10, 10))
fig.suptitle('Network: '+graph[:-4], fontsize=16)

for plotiter in range(1,3):
    if plotiter == 1:
        mu = [0.7]*Nnodes  # 0.8
    else:
        file = open('./DeltasYalmip/delta'+graph[:-4]+'.txt','r')
        mu = []
        for i in range(Nnodes):
            mu.append(float(file.readline()[:-1]))

    ro_over_beta = []
    ro_over_time = np.zeros((1,Tmax))
    itt_beta = 0
    for beta in posible_betas:
        p_Nrep_time = np.zeros((1,Tmax))

        init_array = list(G.nodes)
        random.shuffle(init_array)
        infected = init_array[0:round(ro_init * Nnodes)]
        suscepti = [node for node in init_array if node not in infected]

        for tstep in range(Tmax):
            p_Nrep_time[0,tstep] = len(infected)/Nnodes
            recovered = []
            nowinfect = []

            # 1. For each infected node at time step t, we recover it with probability mu:
            for ninf in infected:
                if random.uniform(0, 1) < mu[int(ninf)-1]:
                    recovered.append(ninf)

            # 2. For each susceptible node at time step t, we traverse all of its neighbors
            # For each infected neighbor (at time step t), the reference node becomes infected with probability beta
            for susnode in suscepti:
                neighbours = list(G.adj[susnode])
                infected_nei = [infnei for infnei in neighbours if infnei in infected]
                counter = 0
                enter = True
                while (enter is True) and (counter < len(infected_nei)):
                    if random.uniform(0, 1) < beta*(G._adj[susnode][infected_nei[counter]]['weight']):
                        nowinfect.append(susnode)
                        enter = False
                    counter = counter + 1

            # 3. Actualize infected() and suscepti()
            for item in recovered:
                infected.remove(item)
            infected = infected + nowinfect
            suscepti = [node for node in init_array if node not in infected]

        ro_over_time[itt_beta, :] = p_Nrep_time[0, :]
        itt_beta = itt_beta + 1
        ro_over_beta.append(np.mean(p_Nrep_time[:,Ttrans:]))

    # Plot one repetition of ro over time for multiple betas
    ax = fig.add_subplot(220+plotiter)
    for iteration in range(ro_over_time.shape[0]):
        ax.plot(ro_over_time[iteration,:], label='Beta(i,j) = '+'w(i,j)*'+str(bet))
    ax.plot([ro_over_beta[0]]*len(ro_over_time[iteration, :]), label='Mean Beta in the stationary')
    ax.legend()
    if np.mod(plotiter,2) != 0:
        ax.set_ylabel('Fraction of infected nodes', fontsize=14)
    ax.set_xlabel('Time t', fontsize=14)
    ax.grid(color='gray', linestyle='-', linewidth=1)
    ax.axis([-10, Tmax, -0.01, 0.3])
    if plotiter == 1:
        ax.title.set_text('Non-optimal curing policy')
    else:
        ax.title.set_text('Optimal curing policy')

    ax = fig.add_subplot(222+plotiter)
    x = np.linspace(1,Nnodes,Nnodes)
    y1 = [0]*Nnodes
    y2 = mu
    ax.plot(x, y2, linewidth=0.3)
    ax.fill_between(x, y1, y2, where=y2>y1, facecolor='#1f77b4', alpha=1)
    ax.plot([np.mean(mu)]*Nnodes, label='Mean Recovery rate')
    ax.legend()
    if plotiter == 1:
        ax.title.set_text('Recovery rate for non-optimal curing policy')
    else:
        ax.title.set_text('Recovery rate for optimal curing policy')
    ax.axis([0, Nnodes, 0, 1])
    if np.mod(plotiter,2) != 0:
        ax.set_ylabel('Recovery rate', fontsize=14)
    ax.set_xlabel('Node', fontsize=14)

fig.savefig('./ImagesSiS/'+'['+graph+']'+'r'+str(ro_init)+'ro_over_time'+'.png')
