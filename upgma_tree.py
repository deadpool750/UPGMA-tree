from Bio.Phylo.BaseTree import Clade, Tree
from Bio import Phylo

def upgma(matrix, labels):
    """
    Performs UPGMA clustering on a given distance matrix.
    Returns a tree dictionary and root label.
    """
    import copy

    current_matrix = copy.deepcopy(matrix)
    current_labels = labels[:]
    tree = {}

    while len(current_labels) > 1:
        # Step 1: Find the closest pair
        min_dist = float("inf")
        i_min, j_min = -1, -1
        for i in range(len(current_matrix)):
            for j in range(i):
                if current_matrix[i][j] < min_dist:
                    min_dist = current_matrix[i][j]
                    i_min, j_min = i, j

        # Step 2: Create new merged label
        label_i = current_labels[i_min]
        label_j = current_labels[j_min]
        new_label = f"({label_j},{label_i})"
        merge_distance = min_dist / 2
        tree[new_label] = (label_j, label_i, merge_distance)

        # Step 3: Compute new row (average distance to others)
        new_row = []
        for k in range(len(current_matrix)):
            if k != i_min and k != j_min:
                d = (current_matrix[i_min][k] + current_matrix[j_min][k]) / 2
                new_row.append(d)

        # Step 4: Build new matrix and label list
        new_matrix = []
        new_labels = []
        for k in range(len(current_matrix)):
            if k not in [i_min, j_min]:
                row = []
                for l in range(len(current_matrix)):
                    if l not in [i_min, j_min]:
                        row.append(current_matrix[k][l])
                new_matrix.append(row)
                new_labels.append(current_labels[k])

        # Add new cluster
        for idx, row in enumerate(new_matrix):
            row.append(new_row[idx])
        new_row.append(0.0)
        new_matrix.append(new_row)
        new_labels.append(new_label)

        current_matrix = new_matrix
        current_labels = new_labels

    if current_labels:
        return tree, current_labels[0]
    else:
        raise ValueError("Tree construction failed: no clusters remaining.")


def build_tree(tree_dict, root_label):
    """
    Builds a Biopython Tree with branch lengths from a UPGMA tree dictionary.
    """
    nodes = {}

    def get_node(label):
        if label in nodes:
            return nodes[label]
        elif label not in tree_dict:
            # Leaf node
            return Clade(branch_length=0, name=label)
        else:
            left_label, right_label, branch_length = tree_dict[label]
            left_node = get_node(left_label)
            right_node = get_node(right_label)
            left_node.branch_length = branch_length
            right_node.branch_length = branch_length
            return Clade(clades=[left_node, right_node])

    return Tree(root=get_node(root_label))


def show_tree(tree):
    """
    Prints the tree as ASCII.
    """
    Phylo.draw_ascii(tree)


def print_branch_lengths(tree):
    """
    Prints branch lengths for all clades with labels.
    """
    print("\nBranch lengths:")
    for clade in tree.find_clades(order="preorder"):
        if clade.name:
            print(f"{clade.name}: branch length = {clade.branch_length}")
