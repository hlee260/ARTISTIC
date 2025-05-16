import re
import os
import pandas as pd
import RNA
import matplotlib.pyplot as plt
from draw_rna.ipynb_draw import draw_struct

# === Helpers ===
def sanitize_sequence(raw: str) -> str:
    return ''.join(re.findall(r"[ATCGU]", raw.upper()))

def highest_salt(buf: str):
    """
    From a buffer‐description string, return (salt_name, concentration_mM).
    Special‐case PBS is converted to to 150 mM NaCl.
    Recognizes NaCl, KCl, MgCl2, and CaCl2 in any order.
    """
    # 0) normalize unicode subscript
    txt = buf.strip().replace("₂", "2")

    # 1) PBS / Phosphate fallback
    if re.search(r"\bPBS\b", txt, re.IGNORECASE) \
    or re.search(r"SALINE", txt, re.IGNORECASE):
        return "NACL", 150.0

    # 2) regex to catch both “Salt … mM” and “… mM Salt”
    pattern = re.compile(r"""
      (?:                                    # either:
        (NaCl|KCl|MgCl2|CaCl2)               #   group1=Salt
        [\s,;]*?                             #   optional sep
        (\d+\.?\d*)\s*mM                     #   group2=conc
      )
      |
      (?:                                    # or:
        (\d+\.?\d*)\s*mM                     #   group3=conc
        [\s,;]*?                             #   optional sep
        (NaCl|KCl|MgCl2|CaCl2)               #   group4=Salt
      )
    """, re.IGNORECASE | re.VERBOSE)

    matches = pattern.findall(txt)
    if not matches:
        return None, None

    # 3) build candidate list
    candidates = []
    for g1, g2, g3, g4 in matches:
        if g1:        # salt first
            salt = g1.upper().replace(" ", "")
            conc = float(g2)
        else:         # conc first
            salt = g4.upper().replace(" ", "")
            conc = float(g3)
            candidates.append((salt, conc))
    # 4) choose highest concentration
    best = max(candidates, key=lambda x: x[1])
    return best

# === Aptamer Database Loader ===
class AptamerDatabase:
    def __init__(self, excel_path: str):
        df = pd.read_excel(excel_path, sheet_name=0)
        df.columns = [
            c.strip().lower().replace(" ", "_")
             .replace("(", "").replace(")", "")
             .replace("/", "_")
            for c in df.columns
        ]
        df['target_norm'] = (
            df['target']
              .astype(str)
              .str.lower()
              .str.replace(r'[^a-z0-9]', '', regex=True)
        )
        self.df = df

    def search_by_target(self, target: str) -> pd.DataFrame:
        # normalize the user’s query the same way
        q = re.sub(r'[^a-z0-9]', '', target.lower())
        return self.df[self.df['target_norm'].str.contains(q)]


# === dART Generator Settings ===
Prom_nt = "TTCTAATACGACTCACTATA"
Prom_t  = "TATAGTGAGTCGTATTAGAA"
O1_nt = "CTACATCCACATACTAATTAAC"
O1_t  = "GTTAATTAGTATGTGGATGTAG"

# Alternative outputs
O2_nt = "CTACTTTCACTTCACAACATCA"
O2_t  = "TGATGTTGTGAAGTGAAAGTAG"
O3_nt = "TACCATCACATTCAATAATCCT"
O3_t  = "AGGATTATTGAATGTGATGGTA"

# Designing Inverter dART for the digital biosensor (O1 output)
O1c_nt = "GTTAATTAGTATGTGGAT"
O1c_t  = "ATCCACATACTAATTAAC"

default_insulation        = "GGGATG"
default_ins_comp          = "CATCCC"
alternative_insulations   = ["GGGAGT","GGGAGA","GGGAAA"]
alternative_ins_comps      = ["ACTCCC","TCTCCC","TTTCCC"]

# Other insulation domains considered if the candidates above do not work!

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
