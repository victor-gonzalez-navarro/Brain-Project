import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import collections
import numpy as np

from sklearn.linear_model import LinearRegression
from networkx.utils import accumulate

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


def plot_degree_distribution(G, title, ax):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    ax.bar(deg, cnt, width=0.01*(max(deg)-min(deg)), color='b')
    ax.title.set_text("Degree Histogram "+title)
    ax.set_ylabel("Count")
    ax.set_xlabel("Degree")


def plot_degree_distribution_line(G,title):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    fig, ax = plt.subplots()
    plt.plot(deg, cnt, color='b')
    plt.title("Degree Histogram "+title)
    plt.ylabel("Count")
    plt.xlabel("Degree")


def plot_adjacencyMatrix(G, title, fig, ax):
    Nnodes = len(G.nodes)
    y, x = np.meshgrid(np.linspace(1, Nnodes, Nnodes), np.linspace(1, Nnodes, Nnodes))
    z = np.array(nx.to_numpy_matrix(G))
    z_min, z_max = z.min(), z.max()
    c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.title.set_text('Adjacency Matrix '+title)
    ax.axis([x.min(), x.max(), y.max(), y.min()])
    fig.colorbar(c, ax=ax)


def plot_binary_adjacencyMatrix(G, title, ax):
    Nnodes = len(G.nodes)
    y, x = np.meshgrid(np.linspace(1, Nnodes, Nnodes), np.linspace(1, Nnodes, Nnodes))
    z = np.array(nx.to_numpy_matrix(G))
    z = z>0
    z_min, z_max = z.min(), z.max()
    ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.title.set_text('Binary Adjacency Matrix '+title)
    ax.axis([x.min(), x.max(), y.max(), y.min()])


def plot_binary_adjacencyMatrix2(G, z, title, fig, ax):
    Nnodes = len(G.nodes)
    y, x = np.meshgrid(np.linspace(1, Nnodes, Nnodes), np.linspace(1, Nnodes, Nnodes))
    z_min, z_max = z.min(), z.max()
    c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.title.set_text('Binary Adjacency Matrix '+title)
    ax.axis([x.min(), x.max(), y.max(), y.min()])
    fig.colorbar(c, ax=ax)


def rich_club_coefficient_vic(G, normalized=True, Q=100):

    if nx.number_of_selfloops(G) > 0:
        raise Exception('rich_club_coefficient is not implemented for '
                        'graphs with self loops.')
    rc = _compute_rc_vic(G)
    if normalized:
        R = G.copy()
        E = R.number_of_edges()
        nx.double_edge_swap(R, Q * E, max_tries=Q * E * 10)
        rcran = _compute_rc_vic(R)
        for it in range(len(rcran)):
            if rcran[it] == 0:
                rcran[it] = 1
        rc = {k: v / rcran[k] for k, v in rc.items()}
    return rc


def _compute_rc_vic(G):
    deghist = nx.degree_histogram(G)
    total = sum(deghist)
    nks = (total - cs for cs in accumulate(deghist) if total - cs > 1)
    edge_degrees = sorted((sorted(map(G.degree, e)) for e in G.edges()),
                          reverse=True)
    ek = G.number_of_edges()
    k1, k2 = edge_degrees.pop()
    rc = {}
    for d, nk in enumerate(nks):
        while k1 <= d:
            if len(edge_degrees) == 0:
                ek = 0
                break
            k1, k2 = edge_degrees.pop()
            ek -= 1
        rc[d] = 2 * ek / (nk * (nk - 1))
    return rc


def compute_correlation(array1, array2, title, ax):
    pcorr = np.corrcoef(array1, array2)[0, 1]
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(array1.reshape(-1, 1), array2)  # perform linear regression
    Y_pred = linear_regressor.predict(array1.reshape(-1, 1))  # make predictions
    ax.plot(array1, array2, 'b.')
    ax.plot(array1, Y_pred, 'r-')
    ax.title.set_text(title +'\nPearson correlation = ' + str(round(pcorr, 3)))


def rich_club_coeff(G, title, ax):
    rich_club_coeff = rich_club_coefficient_vic(G, normalized=True, Q=100)
    vectorx = []
    vectory = []
    for key, val in rich_club_coeff.items():
        vectorx = vectorx + [key, ]
        vectory = vectory + [val, ]
    ax.plot(vectorx, vectory, 'b.')
    ax.title.set_text(title + '\nRich Club coefficient')

# ----------------------------------------------------------------------------------------------------------------- MAIN

# Load Coactivation matrix and Restig-State Connection matrix
Gcoact = nx.read_pajek('./NetworksPajek/Coactivation.net')
Gresti = nx.read_pajek('./NetworksPajek/RestingState.net')
# Both Gcoact and Gresti have the same number of vertices and edges
Grand = nx.gnm_random_graph(len(Gcoact), Gcoact.number_of_edges())


# Plot degree distribution
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131)
plot_degree_distribution(Gcoact, 'Coactivation', ax1)
ax2 = fig.add_subplot(132)
plot_degree_distribution(Gresti, 'Resting-State', ax2)
ax3 = fig.add_subplot(133)
plot_degree_distribution(Grand, 'Random', ax3)


# Plot Adjacency Matrices
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131)
plot_adjacencyMatrix(Gcoact, 'Coactivation', fig, ax1)
ax2 = fig.add_subplot(132)
plot_adjacencyMatrix(Gresti, 'Resting-State', fig, ax2)
ax3 = fig.add_subplot(133)
array3 = np.zeros((len(Grand),len(Grand)))
for item in Grand.edges:
    array3[item[0],item[1]] = 1
    array3[item[1], item[0]] = 1
plot_binary_adjacencyMatrix2(Grand, array3, 'Random', fig, ax3)
plt.show()


# Correlation between Coactivation and Resting-State
array1 = np.array(np.matrix(np.array(nx.to_numpy_matrix(Gcoact))).flatten())[0]
array2 = np.array(np.matrix(np.array(nx.to_numpy_matrix(Gresti))).flatten())[0]
array3 = array3.flatten()
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131)
compute_correlation(array1,array2, 'Coactivation vs Resting-State', ax1)
# Correlation between Coactivation and Random
ax2 = fig.add_subplot(132)
compute_correlation(array1,array3, 'Coactivation vs Random', ax2)
# Correlation between Rsting-State and Random
ax3 = fig.add_subplot(133)
compute_correlation(array2,array3, 'Restin-State vs Random', ax3)
plt.show()


# The rich-club coefficient is a metric on graphs and networks, designed to measure the extent to which well-connected
# nodes also connect to each other. Networks which have a relatively high rich-club coefficient are said to demonstrate
# the rich-club effect and will have many connections between nodes of high degree.
# The rich-club coefficient of a network is useful as a heuristic measurement of the robustness of a network. A high
# rich-club coefficient implies that the hubs are well connected, and global connectivity is resilient to any one hub being removed.
# Paper: https://arxiv.org/pdf/physics/0701290.pdf
# Rich Club Coefficient Coactivation
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131)
rich_club_coeff(nx.Graph(Gcoact), 'Coactivation', ax1)
# Rich Club Coefficient Resting-State
ax2 = fig.add_subplot(132)
rich_club_coeff(nx.Graph(Gresti), 'Restin-State', ax2)
# Rich Club Coefficient Random
ax3 = fig.add_subplot(133)
rich_club_coeff(Grand, 'Random', ax3)
plt.show()






