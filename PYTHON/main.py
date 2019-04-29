import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import collections
import numpy as np



# ------------------------------------------------------------------------------------------------------------ FUNCTIONS
# Ploting coordinates network
def network_plot_3D(G, angle):
    n = G.number_of_nodes()
    deg_vector = list(G.degree)
    deg_list = [deg_vector[i][1] for i in range(n)]
    # Define color range proportional to number of edges adjacent to a single node
    maximumd = max(deg_list)
    colors = [item/maximumd for item in deg_list]

    # 3D network plot
    with plt.style.context(('ggplot')):
        fig = plt.figure(figsize=(10, 7))
        ax = Axes3D(fig)
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for key, value in G._node.items():
            xi = value['x']
            yi = value['y']
            zi = float(value['shape'])
            # Scatter plot
            #ax.scatter(xi, yi, zi, c=colors[int(key)], s=20+20*deg_list[int(key)], edgecolors='k', alpha=0.7)
            ax.scatter(xi, yi, zi, c='r')
            ax.set_xlim([-120, 120])
            ax.set_ylim([-120, 120])
            ax.set_zlim([-120, 120])


    # Set the initial view
    ax.view_init(30, angle)

    plt.show()
    return

def plot_degree_distribution(G):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color='b')
    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticklabels(deg)

def plot_adjacencyMatrix(G):
    Nnodes = len(G.nodes)
    y, x = np.meshgrid(np.linspace(1, Nnodes, Nnodes), np.linspace(1, Nnodes, Nnodes))
    z = np.array(nx.to_numpy_matrix(G))
    z_min, z_max = z.min(), z.max()
    fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.set_title('Adjacency Matrix')
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)


# ----------------------------------------------------------------------------------------------------------------- MAIN

# Load Coactivation matrix and Restig-State Connection matrix
Gcoact = nx.read_pajek('./NetworksPajek/Coactivation.net')
Gresti = nx.read_pajek('./NetworksPajek/RestingState.net')

# Plot degree distribution
plot_degree_distribution(Gcoact)
plot_degree_distribution(Gresti)

# Plot Adjacency Matrices
plot_adjacencyMatrix(Gcoact)
plot_adjacencyMatrix(Gresti)

# Correlation between Coactivation and Resting-State
array1 = np.array(np.matrix(np.array(nx.to_numpy_matrix(Gcoact))).flatten())
array2 = np.array(np.matrix(np.array(nx.to_numpy_matrix(Gresti))).flatten())
pcorr = np.corrcoef(array1, array2)[0, 1]
print('The Pearson Correlation is '+str(pcorr))
#plt.plot(array1, array2, 'ro')

# Display images
plt.show()



