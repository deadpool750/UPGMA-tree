from needleman import pairwise_distance_matrix
from upgma_tree import upgma, build_tree, show_tree, print_branch_lengths
from Bio import SeqIO
from Bio import Phylo

def read_fasta(filename):
    """
    Reads sequences from a FASTA file.

    Parameters:
        filename (str): Path to the FASTA file.

    Returns:
        list[str]: A list of sequences as strings.
    """
    return [str(record.seq) for record in SeqIO.parse(filename, "fasta")]

def save_report(filename, labels, distance_matrix, tree, match=None, mismatch=None, gap=None, sequences=None):
    """
    Saves the full UPGMA tree report to a text file, including input sequences or matrix,
    scoring parameters, ASCII tree view, and branch lengths.

    Parameters:
        filename (str): Output filename.
        labels (list[str]): List of sequence or cluster labels.
        distance_matrix (list[list[float]]): The distance matrix used.
        tree (Tree): Biopython Phylo tree object.
        match (int, optional): Match score (if using sequences).
        mismatch (int, optional): Mismatch penalty (if using sequences).
        gap (int, optional): Gap penalty (if using sequences).
        sequences (list[str], optional): Original sequences, if applicable.
    """
    with open(filename, 'w') as f:
        f.write("UPGMA Tree Report\n")
        f.write("=================\n\n")

        #include sequences if provided, otherwise write matrix
        if sequences:
            f.write("Input Sequences:\n")
            for label, seq in zip(labels, sequences):
                f.write(f"{label}: {seq}\n")
            f.write("\nScoring Parameters:\n")
            f.write(f"  Match: {match}\n")
            f.write(f"  Mismatch: {mismatch}\n")
            f.write(f"  Gap: {gap}\n")
        else:
            f.write("Input Distance Matrix:\n")
            f.write("\t" + "\t".join(labels) + "\n")
            for label, row in zip(labels, distance_matrix):
                f.write(label + "\t" + "\t".join(f"{val:.2f}" for val in row) + "\n")

        # Tree display
        f.write("\nGenerated UPGMA Tree (ASCII View):\n")
        f.write("----------------------------------\n")
        from io import StringIO
        handle = StringIO()
        Phylo.draw_ascii(tree, file=handle)
        f.write(handle.getvalue())

        # Branch lengths
        f.write("\nBranch Lengths:\n")
        for clade in tree.find_clades(order="preorder"):
            if clade.name:
                f.write(f"{clade.name}: branch length = {clade.branch_length}\n")

    print(f"\nReport saved to '{filename}'.")

def main():
    """
    Main entry point for the UPGMA tree builder.
    Allows input of sequences (manually or from FASTA) or a distance matrix.
    Constructs the UPGMA tree, shows it in ASCII, and optionally saves a report.
    """
    print("UPGMA Phylogenetic Tree Builder")

    #select mode of input
    mode = input("Choose input mode:\n1 - Use sequences\n2 - Use distance matrix file\nEnter 1 or 2: ").strip()

    match = mismatch = gap = None
    sequences = None

    #mode 1: Sequences
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

        #get scoring scheme
        match = int(input("Match score: "))
        mismatch = int(input("Mismatch penalty: "))
        gap = int(input("Gap penalty: "))

        #label sequences
        labels = [f"Seq{i+1}" for i in range(len(sequences))]

        #compute distance matrix via Needleman-Wunsch
        distance_matrix = pairwise_distance_matrix(sequences, match, mismatch, gap)

    #mode 2: Distance matrix input
    elif mode == "2":
        matrix_file = input("Enter the distance matrix file (e.g., distance_matrix.txt): ").strip()
        label_file = input("Enter the label file (e.g., labels.txt): ").strip()

        try:
            with open(matrix_file, 'r') as f:
                distance_matrix = [list(map(float, line.strip().split())) for line in f]

            with open(label_file, 'r') as f:
                labels = [line.strip() for line in f]

            #check matrix consistency
            n = len(labels)
            if len(distance_matrix) != n or any(len(row) != n for row in distance_matrix):
                raise ValueError("Matrix dimensions do not match number of labels.")

        except Exception as e:
            print(f"Error reading matrix or labels: {e}")
            return

    else:
        print("Invalid option. Choose 1 or 2.")
        return

    #construct and display UPGMA tree
    try:
        tree_data, root_label = upgma(distance_matrix, labels)
        tree = build_tree(tree_data, root_label)

        print("\nGenerated UPGMA Tree:\n")
        show_tree(tree)
        print_branch_lengths(tree)

        #saving report to text file
        save = input("\nWould you like to save the UPGMA tree and input data to a text report? (y/n): ").strip().lower()
        if save == "y":
            filename = input("Enter filename (e.g., upgma_report.txt): ").strip()
            save_report(
                filename=filename,
                labels=labels,
                distance_matrix=distance_matrix,
                tree=tree,
                match=match,
                mismatch=mismatch,
                gap=gap,
                sequences=sequences
            )

    except Exception as e:
        print(f"Tree generation failed: {e}")

if __name__ == "__main__":
    main()
