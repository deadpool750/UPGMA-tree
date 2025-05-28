from needleman import pairwise_distance_matrix
from upgma_tree import upgma, build_tree, show_tree, print_branch_lengths
from Bio import SeqIO
from Bio import Phylo

def read_fasta(filename):
    return [str(record.seq) for record in SeqIO.parse(filename, "fasta")]

def save_newick(tree, filename="upgma_tree.nwk"):
    Phylo.write(tree, filename, "newick")
    print(f"\nTree saved in Newick format to '{filename}'.")

def main():
    print("UPGMA Phylogenetic Tree Builder")
    mode = input("Choose input mode:\n1 - Use sequences\n2 - Use distance matrix file\nEnter 1 or 2: ").strip()

    if mode == "1":
        strand_source = input("Choose sequence input:\n1 - Load from FASTA file\n2 - Enter manually\nEnter 1 or 2: ").strip()

        if strand_source == "1":
            fasta = input("Enter FASTA file name: ").strip()
            try:
                sequences = read_fasta(fasta)
            except Exception as e:
                print(f"Error reading FASTA: {e}")
                return

        elif strand_source == "2":
            sequences = []
            print("Enter sequences one per line (A, C, G, T, U only). Press Enter on an empty line to finish.")
            while True:
                seq = input(f"Sequence {len(sequences)+1}: ").strip().upper()
                if not seq:
                    break
                if all(c in "ACGTU" for c in seq):
                    sequences.append(seq)
                else:
                    print("Invalid characters detected. Use only A, C, G, T, U.")

        else:
            print("Invalid sequence input method.")
            return

        if len(sequences) < 2:
            print("At least two sequences are required.")
            return

        match = int(input("Match score: "))
        mismatch = int(input("Mismatch penalty: "))
        gap = int(input("Gap penalty: "))

        labels = [f"Seq{i+1}" for i in range(len(sequences))]
        distance_matrix = pairwise_distance_matrix(sequences, match, mismatch, gap)

    elif mode == "2":
        matrix_file = input("Enter the distance matrix file (e.g., distance_matrix.txt): ").strip()
        label_file = input("Enter the label file (e.g., labels.txt): ").strip()

        try:
            with open(matrix_file, 'r') as f:
                distance_matrix = [list(map(float, line.strip().split())) for line in f]

            with open(label_file, 'r') as f:
                labels = [line.strip() for line in f]

            n = len(labels)
            if len(distance_matrix) != n or any(len(row) != n for row in distance_matrix):
                raise ValueError("Matrix dimensions do not match number of labels.")

        except Exception as e:
            print(f"Error reading matrix or labels: {e}")
            return

    else:
        print("Invalid option. Choose 1 or 2.")
        return

    # Build and show tree
    try:
        tree_data, root_label = upgma(distance_matrix, labels)
        tree = build_tree(tree_data, root_label)

        print("\nGenerated UPGMA Tree:\n")
        show_tree(tree)

        print_branch_lengths(tree)

        # Ask to save the tree
        save = input("\nWould you like to save the tree to a Newick file? (y/n): ").strip().lower()
        if save == "y":
            filename = input("Enter filename (e.g., my_tree.nwk): ").strip()
            save_newick(tree, filename)
    except Exception as e:
        print(f"Tree generation failed: {e}")

if __name__ == "__main__":
    main()
