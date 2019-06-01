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


def clustering_coeff_vic(AG, title, ax1, ax2, iterations, color):
    resultclusteringcoef = []
    resultglobaleficienc = []
    numberedges = []
    for i in range(iterations):
        Gaux = nx.from_numpy_matrix(AG)
        resultclusteringcoef.append(nx.average_clustering(Gaux))
        resultglobaleficienc.append(nx.global_efficiency(Gaux))
        numberedges.append(Gaux.number_of_edges())
        minval = np.min(AG[np.nonzero(AG)])
        a, b = np.where(AG == minval)
        if len(a) > 0:
            a = a[0]
        if len(b) > 0:
            b = b[0]
            AG[a, b] = 0
            AG[b, a] = 0
    ax1.plot(numberedges, resultclusteringcoef, color+'-', label=title)
    ax1.set_xlim([min(numberedges), max(numberedges)])
    ax1.title.set_text('Clustering coefficient')
    ax1.set_xlabel('Number of edges')
    ax1.invert_xaxis()
    ax1.grid()
    ax2.plot(numberedges, resultglobaleficienc, color+'-', label=title)
    ax2.set_xlim([min(numberedges), max(numberedges)])
    ax2.title.set_text('Global Efficiency')
    ax2.set_xlabel('Number of edges')
    ax2.invert_xaxis()
    ax2.grid()


# ----------------------------------------------------------------------------------------------------------------- MAIN

# Load Coactivation matrix and Restig-State Connection matrix
Gcoact = nx.read_pajek('./NetworksPajek/Coactivation.net')
Gresti = nx.read_pajek('./NetworksPajek/RestingState.net')
# Both Gcoact and Gresti have the same number of vertices and edges
Grand = nx.read_pajek('./NetworksPajek/Random.net')
#Grand = nx.gnm_random_graph(len(Gcoact), Gcoact.number_of_edges())
#nx.write_pajek(Grand, "./NetworksPajek/Random.net")


# Plot degree distribution
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131)
plot_degree_distribution(Gcoact, 'Coactivation', ax1)
ax2 = fig.add_subplot(132)
plot_degree_distribution(Gresti, 'Resting-State', ax2)
ax3 = fig.add_subplot(133)
plot_degree_distribution(Grand, 'Random', ax3)
fig.savefig('./Images/DegreeDistributions' + '.png')



# Plot Adjacency Matrices
fig = plt.figure(figsize=(16, 4))
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
fig.savefig('./Images/AdjacencyMatrices' + '.png')


# Pearson Correlation between Coactivation and Resting-State
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
fig.savefig('./Images/PearsonCorrelation' + '.png')

# Rich Club Coefficient, Paper: https://arxiv.org/pdf/physics/0701290.pdf
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
fig.savefig('./Images/RichClubCoefficient' + '.png')


# Clustering Coefficient
fig = plt.figure(figsize=(12, 5))
iteration = 5000
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
clustering_coeff_vic(np.array(nx.to_numpy_matrix(Gcoact)), 'Coactivation', ax1, ax2, iteration, 'b')
clustering_coeff_vic(np.array(nx.to_numpy_matrix(Gresti)), 'Resting-State', ax1, ax2, iteration, 'g')
clustering_coeff_vic(array3, 'Random', ax1, ax2, iteration, 'r')
ax1.legend()
fig.savefig('./Images/ClusteringAndEfficiency' + '.png')








# Our brain is a network. It consists of spatially distributed, but functionally linked regions that continuously share
# information with each other. Interestingly, recent advances in the acquisition and analysis of functional neuroimaging
#  data have catalyzed the exploration of functional connectivity in the human brain. Functional connectivity is defined
#  as the temporal dependency of neuronal activation patterns of anatomically separated brain regions and in the past
# years an increasing body of neuroimaging studies has started to explore functional connectivity by measuring the level
#  of co-activation of resting-state fMRI time-series between brain regions. These studies have revealed interesting new
#  findings about the functional connections of specific brain regions and local networks, as well as important new
#  insights in the overall organization of functional communication in the brain network. Here we present an overview of
#  these new methods and discuss how they have led to new insights in core aspects of the human brain, providing an
# overview of these novel imaging techniques and their implication to neuroscience. We discuss the use of spontaneous
# resting-state fMRI in determining functional connectivity, discuss suggested origins of these signals, how functional
# connections tend to be related to structural connections in the brain network and how functional brain communication
#  may form a key role in cognitive performance. Furthermore, we will discuss the upcoming field of examining functional
#  connectivity patterns using graph theory, focusing on the overall organization of the functional brain network.
# Specifically, we will discuss the value of these new functional connectivity tools in examining believed connectivity
# diseases, like Alzheimer's disease, dementia, schizophrenia and multiple sclerosis.


