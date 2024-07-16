import matplotlib.pyplot as plt
import networkx as nx


# Function to create and save the graph
def create_network_graph(df, betti_number):
    # Filter DataFrame for the specified Betti number
    filtered_df = df[df["betti_number"] == betti_number]

    # Create a NetworkX graph
    G = nx.Graph()

    # Add nodes and edges to the graph
    for gene_list in filtered_df["gene_set"]:
        for i in range(len(gene_list)):
            for j in range(i + 1, len(gene_list)):
                gene1, gene2 = gene_list[i], gene_list[j]
                if G.has_edge(gene1, gene2):
                    G[gene1][gene2]["weight"] += 0.5
                else:
                    G.add_edge(gene1, gene2, weight=1)

    # Function to save the graph
    def save_graph(G, title="Gene Interaction Network (Betti {})".format(betti_number)):
        plt.figure(figsize=(12, 12))
        pos = nx.spring_layout(G)

        # Get edge weights for width
        edges = G.edges(data=True)
        edge_widths = [data["weight"] for _, _, data in edges]

        # Draw the graph with edge widths
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=50,
            font_size=10,
            font_weight="bold",
            width=edge_widths,
            edge_color="gray",
        )
        plt.title(title)
        plt.savefig(f"output/network_graphs/{title}.png")
        plt.close()

    # Save the graph
    save_graph(G)
