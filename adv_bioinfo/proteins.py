import time
from Bio import Entrez
import pandas as pd

# Identify yourself to NCBI
Entrez.email = "d00447258@utahtech.edu"

# Load accession list (from CSV or acc_list.txt)
accs = pd.read_csv("CovidAccessionNumbers.csv")["acc_numbers"].tolist()

with open("covid_proteins.fasta", "w") as out_f:
    for i, acc in enumerate(accs, 1):
        try:
            # Map nucleotide accession -> protein records
            link = Entrez.elink(dbfrom="nuccore", db="protein", id=acc, linkname="nuccore_protein")
            link_result = Entrez.read(link)

            protein_ids = []
            for l in link_result[0].get("LinkSetDb", []):
                protein_ids.extend([x["Id"] for x in l["Link"]])

            if not protein_ids:
                print(f"[{i}/{len(accs)}] No proteins found for {acc}")
                continue

            # Fetch protein sequences in FASTA
            handle = Entrez.efetch(db="protein", id=",".join(protein_ids),
                                   rettype="fasta", retmode="text")
            seqs = handle.read()
            out_f.write(seqs)

            print(f"[{i}/{len(accs)}] Downloaded {len(protein_ids)} proteins for {acc}")

            # Sleep to avoid NCBI rate-limits
            time.sleep(0.5)

        except Exception as e:
            print(f"[{i}/{len(accs)}] Error with {acc}: {e}")
            time.sleep(2)  # wait longer before retry

