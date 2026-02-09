import pandas as pd
import re

# write into a class
#GFF_PATH = "PROKKA_SYNBICT/PROKKA_SYNBICT.gff"
#BLAST_PATH = "PROKKA_SYNBICT/PROKKA_SYNBICT.proteins.tmp.*.blast"

class ProkkaParser:
    def __init__(self, gff_path, blast_path):
        self.gff_path = gff_path
        self.blast_path = blast_path
        self.merged_df = None

    def parse_gff_and_blast(self):
        lines = []
        with open(self.gff_path) as f:
            for line in f:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                lines.append(line)

        df = pd.read_csv(
            pd.io.common.StringIO("".join(lines)),
            sep="\t",
            header=None,
            names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
        )
        # pull common fields from the attributes string
        df["gene"] = df["attributes"].apply(lambda x: self.get_attr(x, "gene"))
        df["product"] = df["attributes"].apply(lambda x: self.get_attr(x, "product"))
        df["locus_tag"] = df["attributes"].apply(lambda x: self.get_attr(x, "locus_tag"))
        df["ID"] = df["attributes"].apply(lambda x: self.get_attr(x, "ID"))
        ab = df["attributes"].apply(lambda x: self.get_attr(x, "database_protein.fasta"))
        df["protein_name"] = [x.split(":")[4] if x is not None else "" for x in ab]
        df["custom_db"] = [x.split(":")[3] if x is not None else "" for x in ab]
        # before map to blast, first remove the not CDS type
        # if type == CDS
        only_cds_df = df[df["type"] == "CDS"]
        not_cds_df = df[df["type"] != "CDS"]
        # reset index after filtering
        only_cds_df = only_cds_df.reset_index(drop=True)
        only_cds_df["query_index"] = only_cds_df.index
        # Common reasons:
        # Non-CDS features removed (RNA not in protein FASTA)
        # Short CDS filtered
        # Order preserved, but IDs skipped
        # BLAST only reports queries with hits
        
        df_blast = self.parse_blast_to_dataframe()
        # Perform an outer merge on the 'ID' column
        merged_df = only_cds_df.merge(
            df_blast,
            on="query_index",
            how="left"
        )
        
        # non-CDS rows must have ALL merged_df columns
        for col in merged_df.columns:
            if col not in not_cds_df.columns:
                not_cds_df[col] = pd.NA

        # query_index should be NaN for non-CDS
        not_cds_df["query_index"] = pd.NA
        # reorder columns to match merged_df exactly
        not_cds_df = not_cds_df[merged_df.columns]

        # concatenate
        final_df = pd.concat(
            [merged_df, not_cds_df],
            ignore_index=True
        )

        return final_df

    def get_attr(self, attr_str, key):
        if pd.isna(attr_str):
            return None
        for item in str(attr_str).split(";"):
            if item.startswith(key + "=") or item.__contains__(key):
                return item.split("=", 1)[1]
        return None

    def parse_blast_to_dataframe(self):
        """
        Parse BLASTP text output into a pandas DataFrame.

        Columns:
        1) query_index : int   (Query number)
        2) identity_pct: int   (percentage identity)
        3) protein_id      : str   (FASTA header)
        """
        rows = []

        with open(self.blast_path, "r") as f:
            lines = f.readlines()

        current_query = None
        current_header = None
        collecting_header = False

        for line in lines:
            line = line.rstrip()

            # Capture Query=
            if line.startswith("Query="):
                q = line.split("=", 1)[1].strip()
                if q.isdigit():
                    current_query = int(q) - 1 
                else:
                    current_query = None
                current_header = None
                collecting_header = False
                continue

            # Start of FASTA header
            if line.startswith(">"):
                current_header = line[1:].strip()
                collecting_header = True
                continue

            # Wrapped FASTA header lines
            if collecting_header:
                if line.startswith("Length="):
                    collecting_header = False
                elif line.strip():
                    current_header += " " + line.strip()
                continue

            # Identity percentage
            if (
                "Identities =" in line
                and current_header is not None
                and current_query is not None
            ):
                match = re.search(r"\((\d+)%\)", line)
                if match:
                    rows.append(
                        {
                            "query_index": current_query,
                            "identity_pct": int(match.group(1)),
                            "protein_id": current_header.split(" ")[0],
                        }
                    )
                    current_header = None  # reset for next hit

        return pd.DataFrame(rows, columns=["query_index", "identity_pct", "protein_id"])