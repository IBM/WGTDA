import dcor
import numpy as np
import gudhi as gd
from scipy.stats import pearsonr
from tqdm import tqdm
import math
from gudhi.representations import Landscape
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from helper import suppress_stdout_stderr

def set_plot_fonts():
    """
    Sets consistent Helvetica fonts with sizes ≥10 points.
    Adjust these values as needed to meet your thesis formatting guidelines.
    """
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica']  # If not installed, falls back to another sans-serif
    plt.rcParams['axes.titlesize'] = 200
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['legend.title_fontsize'] = 14


def save_persistence_diagram(diag, filename="persistence_diagram.png", title="Persistent Diagram (Cancer Type)", dpi=600):
    
    
    # Plot and customize the persistence diagram
    plt.figure(figsize=(8, 6))
    gd.plot_persistence_diagram(diag, legend=True)
    plt.title(title)
    plt.xlabel("Birth")
    plt.ylabel("Death")
    #plt.grid(alpha=0.5)
    
    # Save the plot as an image
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Persistence diagram saved as {filename}")

def compute_distance_correlation_matrix(DF):
    
    num_genes = DF.shape[1]
    dist = np.zeros((num_genes, num_genes))
    
    for i in range(num_genes):
        
        for j in range(i+1, num_genes):
            
            dist[i,j] = dcor.distance_correlation(DF[:,i], DF[:,j]) #Distance Correlations 
    
    dist = dist + dist.T + np.eye(num_genes)
    
    return  dist

def compute_pearson_correlation_matrix(DF: np.ndarray) -> np.ndarray:
    """
    Compute the Pearson correlation matrix for the given gene expression data.

    Parameters:
    - DF: np.ndarray, gene expression data with rows as samples and columns as genes.

    Returns:
    - corr_matrix: np.ndarray, the Pearson correlation matrix.
    """
    num_genes = DF.shape[1]
    corr_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(i, num_genes):
            if i == j:
                corr_matrix[i, j] = 1.0  # Perfect correlation with itself
            else:
                corr, _ = pearsonr(DF[:, i], DF[:, j])
                corr_matrix[i, j] = corr
                corr_matrix[j, i] = corr  # Symmetric matrix


    return 1- corr_matrix

def patient_correlation_measure(F, M):
    
    F = F.T
    num_genes = M.shape[1]
    dist = np.zeros((num_genes, num_genes))
        
    for i in range(num_genes):
        for j in range(i+1, num_genes):
            
            dist[i,j] = M[i,j] + math.sqrt(F[i] + F[j])/10 #smaller values are shrunk and larger values are inflated 
    
    dist = dist + dist.T + np.eye(num_genes)
    
    return dist

def compute_wto_matrix(gene_exp_arr: np.ndarray, adj_matrix="dc")-> np.ndarray:
    """
    Compute the Weighted Topological Overlap (wTO) matrix.

    Parameters:
    - df: np.ndarray, input data.

    Returns:
    - wto_matrix: np.ndarray, Topological Overlap matrix
    """
    if adj_matrix == "dc":
        adjacency_matrix = compute_distance_correlation_matrix(DF=gene_exp_arr)
    elif adj_matrix == "pearsons":
        adjacency_matrix = compute_pearson_correlation_matrix(DF=gene_exp_arr)
    else:
        raise ValueError("Invalid adjacency matrix type. Choose 'dc' or 'pearsons'.")

    num_genes = adjacency_matrix.shape[0]
    wto_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(num_genes):
            if i != j:
                min_ki_kj = (
                    min(
                        np.sum(np.abs(adjacency_matrix[i, :])),
                        np.sum(np.abs(adjacency_matrix[j, :])),
                    )
                    + 1
                    - np.abs(adjacency_matrix[i, j])
                )
                wto_matrix[i, j] = (
                    np.dot(adjacency_matrix[i, :], adjacency_matrix[:, j])
                    + adjacency_matrix[i, j]
                ) / min_ki_kj

    return wto_matrix

def wgtda(preprocessing: np.ndarray, dimensions=2, gene_dict=None, file_name="", title=""):

    interactions_list = []

    interactions = pd.DataFrame(
        columns=['interaction_id', 'betti_number', 'birth', 'death', 'lifespan', 'birth_nodes',
                 'death_nodes', 'birth_geneset',
                 'death_geneset', 'geneset'])

    """
    Construct a Vietoris-Rips complex from the specified preprocessing matrix and compute its persistent homology.

    Parameters:
    - preprocessing : ndarray
        A square matrix where element [i, j] represents the distance between the i-th and j-th elements.
    - dimension : int, default = 3
        The maximum dimension of simplices to be considered in the Vietoris-Rips complex. Default is 3.

    Returns:
    - tuple[PersistentHomologyComputer, FilteredSimplicialComplex]
        A tuple containing the persistent homology computer and the Vietoris-Rips complex.
        The PersistentHomologyComputer object has methods to access and manipulate the computed homology.
        The FilteredSimplicialComplex object represents the simplicial complex constructed from the given distance matrix.

    Notes:
    - The function uses the 'matilda' library, which must be installed and properly configured in the environment.
    - This function is intended for use in computational topology, particularly in the analysis of high-dimensional data.

    Examples:
    - Constructing the complex and computing homology for a given distance matrix:
         matrix = np.array([[0, 1, 1.5], [1, 0, 2], [1.5, 2, 0]])
     persistence, rips_complex = VRcomplex(matrix, dimension=2)
    """
    # Create a FilteredSimplicialComplex object

    rips_complex = gd.RipsComplex(distance_matrix=preprocessing, max_edge_length=float('Inf')).create_simplex_tree(
                max_dimension=dimensions)  # Weights used include per-patient gene expressions
    rips_complex.collapse_edges()
    rips_complex.expansion(3)
    
    # Construct the Vietoris-Rips complex from the preprocessed distance matrix
    betti_pairs = rips_complex.persistence(persistence_dim_max=dimensions)

    save_persistence_diagram(betti_pairs,filename=file_name, title=title)

    betti_pairs.reverse()
    persistence_pairs = rips_complex.persistence_pairs()



    for index, ((betti_number, (birth, death)), (birth_nodes, death_nodes)) in enumerate(
            zip(betti_pairs, persistence_pairs)):
            birth_genes = [gene_dict[item] for item in birth_nodes]
            death_genes = [gene_dict[item] for item in death_nodes]
            gene_set = birth_genes + death_genes
            lifespan = death - birth
            interaction = [
                index, betti_number, birth, death, lifespan, birth_nodes, death_nodes, birth_genes, death_genes,
                gene_set
            ]
            interactions_list.append(interaction)

             # Create a DataFrame from the NumPy array
    new_interactions = pd.DataFrame(interactions_list, columns=interactions.columns)

    # Concatenate the new interactions with the existing DataFrame
    interactions = pd.concat([interactions, new_interactions], ignore_index=True)

    return interactions

def compute_simplicial_complex_and_landscapes(patient_data: np.ndarray, dist_corr_matrix: np.ndarray, num_landscape=2, resolution=10, dim=2):
    """
    Compute persistent homology and persistence landscapes for all patients.
    
    Parameters:
    - patient_data: np.ndarray, gene expression data with rows as patients and columns as genes.
    - dist_corr_matrix: np.ndarray, population-wide distance correlation matrix.
    - max_filtration: float, maximum filtration value to limit the persistence computation.
    
    Returns:
    - persistence_landscapes: np.ndarray, persistence landscapes for all patients.
    """
    Persistent_diagrams0, Persistent_diagrams1, Persistent_diagrams2 = [], [], []
    #print("Computing the Simplicial Complex and persistence for patients")
    for patient_exp in tqdm(patient_data):
        patient_matrix = patient_correlation_measure(patient_exp, dist_corr_matrix)
        #print(patient_matrix.shape)
        rips_complex = gd.RipsComplex(distance_matrix=patient_matrix, max_edge_length=float('Inf')).create_simplex_tree(
                max_dimension=dim)  # Weights used include per-patient gene expressions
        rips_complex.collapse_edges()
        rips_complex.expansion(3)
        diag = rips_complex.persistence()

        Persistent_diagrams0.append(rips_complex.persistence_intervals_in_dimension(0))
        Persistent_diagrams1.append(rips_complex.persistence_intervals_in_dimension(1))
        Persistent_diagrams2.append(rips_complex.persistence_intervals_in_dimension(2))


        #remove_infinity = lambda barcode: np.array([bars for bars in barcode if bars[1] != np.inf])

    #Persistent_diagrams0 = list(map(remove_infinity, Persistent_diagrams0))
    with suppress_stdout_stderr():
        landscape = Landscape(num_landscapes=num_landscape, resolution=resolution)
        samplelandscape0lnd = landscape.fit_transform(Persistent_diagrams0)
        samplelandscape1lnd = landscape.fit_transform(Persistent_diagrams1)
        samplelandscape2lnd = landscape.fit_transform(Persistent_diagrams2)
        landscapes = np.column_stack((samplelandscape0lnd, samplelandscape1lnd, samplelandscape2lnd))

    return landscapes

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"


