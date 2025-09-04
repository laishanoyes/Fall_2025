from Bio import Entrez
import pandas as pd
import time

# Always tell NCBI who you are
Entrez.email = "your.email@domain.com"

# Load accession numbers from your CSV
accs = pd.read_csv("CovidAccessionNumbers.csv")["acc_numbers"].tolist()

with open("covid_genomes.fasta", "w") as out_f:
    for i, acc in enumerate(accs, 1):
        try:
            # Fetch nucleotide sequence directly
            handle = Entrez.efetch(db="nucleotide", id=acc,
                                   rettype="fasta", retmode="text")
            seq_record = handle.read()
            out_f.write(seq_record)

            print(f"[{i}/{len(accs)}] Downloaded genome for {acc}")

            # Sleep to avoid overloading NCBI servers
            time.sleep(0.5)

        except Exception as e:
            print(f"[{i}/{len(accs)}] Error with {acc}: {e}")
            time.sleep(2)

