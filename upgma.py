import copy

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

def upgma(matrix, labels):
    tree = {}
    current_matrix = copy.deepcopy(matrix)
    current_labels = labels[:]

    while len(current_labels) > 1:
        #1: find closest pair
        min_dist = float("inf")
        i_min, j_min = -1, -1
        for i in range(len(current_matrix)):
            for j in range(i):
                if current_matrix[i][j] < min_dist:
                    min_dist = current_matrix[i][j]
                    i_min, j_min = i, j

        #2: merging labels
        label_i = current_labels[i_min]
        label_j = current_labels[j_min]
        new_label = f"({label_j},{label_i})"  # using j,i to match lower triangle

        # Store tree step
        tree[new_label] = (label_j, label_i, min_dist / 2)

        #3: Computing new row
        new_row = []
        for k in range(len(current_matrix)):
            if k != i_min and k != j_min:
                d = (current_matrix[i_min][k] + current_matrix[j_min][k]) / 2
                new_row.append(d)

        # Step 4: Build new matrix
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

        # Append new row/column for new cluster
        for idx, row in enumerate(new_matrix):
            row.append(new_row[idx])
        new_row.append(0.0)  # distance to self
        new_matrix.append(new_row)
        new_labels.append(new_label)

        # Update matrix and labels
        current_matrix = new_matrix
        current_labels = new_labels

    return tree, current_labels[0]


def build_tree(tree_dict, root_label):
    nodes = {}

    def get_node(label):
        if label in nodes:
            return nodes[label]
        elif label not in tree_dict:
            # It's a leaf node
            node = Clade(branch_length=0, name=label)
        else:
            left, right, length = tree_dict[label]
            left_node = get_node(left)
            right_node = get_node(right)
            node = Clade(branch_length=0, clades=[
                Clade(branch_length=length, name=None, clades=[left_node]),
                Clade(branch_length=length, name=None, clades=[right_node])
            ])
        nodes[label] = node
        return node

    return Tree(root=get_node(root_label))

names = ["A", "B", "C", "D", "E"]
matrix = [
    [0, 16, 6, 16, 6],
    [16, 0, 16, 8, 16],
    [6, 16, 0, 16, 2],
    [16, 8, 16, 0, 16],
    [6, 16, 2, 16, 0]
]

tree_data, root_label = upgma(matrix, names)
bio_tree = build_tree(tree_data, root_label)

Phylo.draw_ascii(bio_tree)



