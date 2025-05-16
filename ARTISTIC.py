import re
import os
import pandas as pd
import RNA
import matplotlib.pyplot as plt
from draw_rna.ipynb_draw import draw_struct

# === Helpers ===
def sanitize_sequence(raw: str) -> str:
    """Keep only valid bases A, T, C, G, U (uppercase)."""
    return ''.join(re.findall(r"[ATCGU]", raw.upper()))

def highest_salt(buf: str):
    """
    From a buffer‐description string, return (salt_name, concentration_mM).
    Special‐case PBS to 150 mM NaCl.
    """
    txt = buf.strip()
    # 1) Special case PBS
    if "PBS" in txt.upper() or "PHOSPHATE BUFFERED SALINE" in txt.upper():
        return "NACL", 150

    # 2) Case‐insensitive match for NaCl, KCl, MgCl2 with concentration in mM
    matches = re.findall(
        r"(NaCl|KCl|MgCl2)\s*,?\s*(\d+\.?\d*)\s*mM",
        txt,
        flags=re.IGNORECASE
    )
    if not matches:
        return None, None

    # 3) Pick the highest concentration
    salt, conc = max(matches, key=lambda x: float(x[1]))
    return salt.upper(), float(conc)


# === Aptamer Database Loader ===
class AptamerDatabase:
    def __init__(self, excel_path: str):
        df = pd.read_excel(excel_path, sheet_name=0)
        df.columns = [
            c.strip().lower()
             .replace(" ", "_")
             .replace("(", "")
             .replace(")", "")
             .replace("/", "_")
            for c in df.columns
        ]
        self.df = df

    def search_by_target(self, target: str) -> pd.DataFrame:
        return self.df[self.df['target'].str.contains(target, case=False, na=False)]


# === dART Generator Settings ===
Prom_nt = "TTCTAATACGACTCACTATA"
Prom_t  = "TATAGTGAGTCGTATTAGAA"
O1_nt = "CTACATCCACATACTAATTAAC"
O1_t  = "GTTAATTAGTATGTGGATGTAG"

# Alternative outputs if O1 output is not compatible
O2_nt = "CTACTTTCACTTCACAACATCA"
O2_t  = "TGATGTTGTGAAGTGAAAGTAG"
O3_nt = "TACCATCACATTCAATAATCCT"
O3_t  = "AGGATTATTGAATGTGATGGTA"

# Designing Inverter dART for the digital biosensor
O1c_nt = "GTTAATTAGTATGTGGAT"
O1c_t  = "ATCCACATACTAATTAAC"

default_insulation        = "GGGATG"
default_ins_comp          = "CATCCC"
alternative_insulations   = ["GGGAGT","GGGAGA","GGGAAA"]
alternative_ins_comps      = ["ACTCCC","TCTCCC","TTTCCC"]
alternative_alt_insulations = [
    "GGGAAC","GGGAAG","GGGAAT","GGGACA","GGGACG","GGGACT","GGGAGC","GGGAGG",
    "GGGATA","GGGATC","GGGATT","GGGCAA","GGGCAC","GGGCAG","GGGCAT","GGGCCA",
    "GGGCCC","GGGCCG","GGGCCT","GGGCGA","GGGCGC","GGGCGG","GGGCGT","GGGCTA",
    "GGGCTC","GGGCTG","GGGCTT","GGGTAA","GGGTAC","GGGTAG","GGGTAT","GGGTCA",
    "GGGTCC","GGGTCG","GGGTCT","GGGTGA","GGGTGC","GGGTGG","GGGTGT","GGGTTA",
    "GGGTTC","GGGTTG","GGGTTT"
]
alternative_alt_comps = [
    "GTTCCC","CTTCCC","ATTCCC","TGTCCC","CGTCCC","TGACCC","GCTCCC","CCTCCC",
    "TATCCC","GATCCC","AATCCC","TTGCCC","GTGCCC","CTGCCC","ATGCCC","TGGCCC",
    "GGGCCC","CGGCCC","AGGCCC","TCGCCC","GCGCCC","CCGCCC","ACGCCC","TAGCCC",
    "GAGCCC","CAGCCC","AAGCCC","TTACCC","GTACCC","CTACCC","ATACCC","TGACCC",
    "GGACCC","CGACCC","AGACCC","TCACCC","GCACCC","CCACCC","ACACCC","TAACCC",
    "GAACCC","CAACCC","AAACCC"
]

Output_nt = O1_nt
Output_t  = O1_t

bp_threshold = 0.9  # base‐pair probability threshold for binding region

def create_dART_template(aptamer: str, insulation: str, insulation_comp: str):
    """Build the DNA template strand for dART."""
    encoded = f"{Output_t}{insulation}{aptamer}{insulation_comp}"
    return f"{encoded}{Prom_t}", insulation, insulation_comp

def get_rna_transcript(encoded_seq: str) -> str:
    """Convert DNA‐encoded domain to its RNA transcript."""
    comp = {'A':'U','T':'A','U':'A','C':'G','G':'C'}
    return ''.join(comp[b] for b in reversed(encoded_seq))

def test_insulations(aptamer_seq: str):
    """
    Test all insulation domains for a given aptamer sequence.
    Returns (True, dART_template, insulation, insulation_comp) on success,
    or (False, None, None, None) if no valid design found.
    """
    apt = sanitize_sequence(aptamer_seq)
    all_pairs = [(default_insulation, default_ins_comp)]
    all_pairs += list(zip(alternative_insulations, alternative_ins_comps))
    all_pairs += list(zip(alternative_alt_insulations, alternative_alt_comps))

    for ins, ins_comp in all_pairs:
        dna_tmpl, ins, ins_c = create_dART_template(apt, ins, ins_comp)
        rna = get_rna_transcript(f"{Output_t}{ins}{apt}{ins_comp}")

        RNA.cvar.temperature = 37
        struct, mfe = RNA.fold(rna)
        fc = RNA.fold_compound(rna)
        fc.pf()
        bpp = fc.bpp()

        start  = range(7)
        target = range(6 + len(apt), 6 + len(apt) + 8)
        count  = sum(1 for i in start for j in target if bpp[i][j] > bp_threshold)

        if count == 6:
            return True, dna_tmpl, ins, ins_comp

    return False, None, None, None
