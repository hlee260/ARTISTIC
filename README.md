# ARTISTIC
ARTISTIC was run with python 3.12.7, but probaboly works with newer versions. Most processes are using the following packages:

import re
import os
import pandas as pd
import RNA
import matplotlib.pyplot as plt
from draw_rna.ipynb_draw import draw_struct

ARTISTIC contains the functions for generating dART sequences and extracting buffer conditions for titration experiment.
Main Script is used to input the target ligand. ARTISTIC find aptamers from UTexas_Datasheet, extract the aptamer name, sequence, buffer conditions, and reported Kd.
ARTISTIC will produce an excel file containing the dART sequences to order, and prescribe experiments for titrating different salts and ligand concentrations.
