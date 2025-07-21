import pandas as pd
import seaborn as sns
import ast
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rcParams

def set_plot_fonts():
    """
    Sets consistent Helvetica fonts with sizes ≥10 points.
    Adjust these values as needed to meet your thesis formatting guidelines.
    """
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica']  # If not installed, falls back to another sans-serif
    rcParams['axes.titlesize'] = 18
    rcParams['axes.labelsize'] = 14
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12
    rcParams['legend.fontsize'] = 12
    rcParams['legend.title_fontsize'] = 14

def plot_gene_network(df, geneset_col='geneset', layout_func=nx.spring_layout, scale_factor=3.0, 
                      node_size=800, figsize=(12,8), title="Gene Interaction Network", save_path=None, color_number=0):
    """
    Plots a gene interaction network from a DataFrame and optionally saves it.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the interactions. Must have a column with lists of genes.
    geneset_col : str
        The column name in df that contains the gene lists.
    layout_func : callable
        The layout function for positioning nodes (e.g., nx.spring_layout).
    scale_factor : float
        A scaling factor for edge widths. 
    node_size : int
        Size of the nodes in the plot.
    figsize : tuple
        Figure size.
    title : str
        Title of the plot.
    save_path : str or None
        Path to save the figure. If None, the plot will not be saved.
    """

    import matplotlib.pyplot as plt
    import networkx as nx
    import ast

    set2_color = sns.color_palette("pastel")[color_number]

    # Convert the geneset column from string to list if needed
    df[geneset_col] = df[geneset_col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Initialize a graph
    G = nx.Graph()

    # Dictionary for edge weights
    edge_weights = {}

    # Build the edges from the DataFrame
    for idx, row in df.iterrows():
        genes = row[geneset_col]
        for i in range(len(genes)):
            for j in range(i+1, len(genes)):
                if genes[i] != genes[j]:
                    edge = tuple(sorted([genes[i], genes[j]]))
                    edge_weights[edge] = edge_weights.get(edge, 0) + 1

    # Add edges to the graph
    for edge, weight in edge_weights.items():
        G.add_edge(*edge, weight=weight)

    # Plot the graph
    plt.figure(figsize=figsize)
    pos = layout_func(G, k=2.0, seed=42)

    edges = G.edges(data=True)
    weights = [attr['weight'] for _, _, attr in edges]
    max_weight = max(weights) if weights else 1
    edge_widths = [scale_factor * (w / max_weight) for w in weights]

    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=set2_color, edgecolors='black')
    nx.draw_networkx_labels(G, pos, font_size=7, font_color="black")
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color="gray")

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()

    # Save the plot if a path is specified
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"Plot saved successfully to {save_path} at 600 DPI.")

    plt.show()

import pandas as pd
import ast
import networkx as nx
from pyvis.network import Network

def plot_gene_network_interactive(df, geneset_col='geneset', output_file='../../index.html'):
    """
    Creates an interactive gene network visualization using pyvis.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the interactions. Must have a column with lists of genes.
    geneset_col : str
        The column name in df that contains the gene lists.
    output_file : str
        The name of the HTML file to store the network visualization.
    """

    # Convert the geneset column from string to list if needed
    df[geneset_col] = df[geneset_col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Initialize a graph
    G = nx.Graph()

    # Dictionary for edge weights
    edge_weights = {}

    # Build the edges from the DataFrame
    for idx, row in df.iterrows():
        genes = row[geneset_col]
        # Create all pairwise combinations
        for i in range(len(genes)):
            for j in range(i+1, len(genes)):
                if genes[i] != genes[j]:
                    edge = tuple(sorted([genes[i], genes[j]]))
                    edge_weights[edge] = edge_weights.get(edge, 0) + 1

    # Add edges to the graph
    for edge, weight in edge_weights.items():
        gene_a, gene_b = edge
        G.add_edge(gene_a, gene_b, weight=weight)

    # Create a pyvis network
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black", notebook=True, cdn_resources="in_line")

    # Convert the networkx graph to pyvis network
    net.from_nx(G)
    # You can customize the physics of the network for better layout
    net.toggle_physics(True)
    
    # Show labels or add hover information
    # By default, hovering over a node shows the node label.
    # You can also customize tooltips or add more info.

    # Save and display
    net.show(output_file)

# Example usage:
# df = pd.read_csv("interactions.csv")
# plot_gene_network_interactive(df, geneset_col='geneset', output_file='network.html')

def plot_gene_network_interactive2(df, geneset_col='geneset', output_file='network.html', title=""):
    """
    Creates an enhanced interactive gene network visualization using pyvis.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the interactions. Must have a column with lists of genes.
    geneset_col : str
        The column name in df that contains the gene lists.
    output_file : str
        The name of the HTML file to store the network visualization.
    """
    # Convert the geneset column from string to list if needed
    df[geneset_col] = df[geneset_col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Initialize a graph
    G = nx.Graph()

    # Dictionary for edge weights
    edge_weights = {}

    # Build the edges from the DataFrame
    for idx, row in df.iterrows():
        genes = row[geneset_col]
        # Create all pairwise combinations
        for i in range(len(genes)):
            for j in range(i+1, len(genes)):
                if genes[i] != genes[j]:
                    edge = tuple(sorted([genes[i], genes[j]]))
                    edge_weights[edge] = edge_weights.get(edge, 0) + 1

    # Add edges and calculate node degrees
    for edge, weight in edge_weights.items():
        gene_a, gene_b = edge
        G.add_edge(gene_a, gene_b, weight=weight, title=f"Weight: {weight}")

    # Compute node degrees for scaling
    node_degrees = dict(G.degree())
    max_degree = max(node_degrees.values()) if node_degrees else 1

    # Add nodes with hover info and scaled sizes
    for node in G.nodes():
        G.nodes[node]['size'] = 10 + (node_degrees[node] / max_degree) * 20  # Scale node size
        G.nodes[node]['title'] = f"Gene: {node}\nConnections: {node_degrees[node]}"
        G.nodes[node]['color'] = 'lightblue'  # Default color

    # Create a pyvis network
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black", notebook=True, cdn_resources="in_line")
    net.from_nx(G)
    net.toggle_physics(True)


    # Generate the HTML and customize it with a heading
    temp_file = "temp_file.html"
    net.write_html(temp_file)

    # Add a heading to the generated HTML
    with open(temp_file, "r") as f:
        html_content = f.read()

    heading_html = f"""
    <h1 style="text-align: center; font-family: Arial, sans-serif; margin-top: 20px;">{title}</h1>
    """
    # Insert the heading before the network div
    html_content = html_content.replace("<body>", f"<body>{heading_html}")

    # Write back the final HTML
    with open(output_file, "w") as f:
        f.write(html_content)

    print(f"Interactive network saved to: {output_file}")

    # Save and display
    net.show(output_file)
