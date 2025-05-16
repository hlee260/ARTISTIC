# ARTISTIC
ARTISTIC was run with python 3.12.7. All processes are using the following packages (please download one-by-one if necessary):

import re

import os

import pandas as pd

import RNA

import matplotlib.pyplot as plt

from draw_rna.ipynb_draw import draw_struct

# Running the code

ARTISTIC contains the functions for generating dART sequences and extracting buffer conditions for titration experiments. Download ARTISTIC.py and UTexas_Datasheet.xlsx [Askari, A., et al., Nuc. Acid Res. 52, 351–359 (2024)] to use them with the main script.

Main was run as a Jupyter Notebook file. Running the code will ask you to input the target ligand. ARTISTIC find aptamers from UTexas_Datasheet, extracts the aptamer name, sequence, buffer conditions, and reported Kd. The example shows building dARTs for thrombin.

ARTISTIC will produce an excel file containing the dART sequences to order, and prescribe experiments for titrating different salts and ligand concentrations.
