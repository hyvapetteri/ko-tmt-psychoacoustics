#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:21:44 2024

@author: hyvario3
"""

import pandas as pd
from itertools import product
from scipy.io import loadmat
from glob import glob

subj_ids = ['10', '13', '14', '24', '25', '28', '29', '19', '26', '31']
attrs = ['Loudness','Reverberance']
reps = [0,1,2,3]

for subj, attr, rep in product(subj_ids, attrs, reps):
    df = pd.read_csv(f"../data/{subj}_{rep}-{attr}.csv")
    df_clean = df.query("COMPLETED == 1")
    df_clean["ANS"] = df_clean["ANS"].astype(int)
    
    df_clean[["SUBID", "ATTR", "KO", "A", "B", "ANS"]].to_csv(f"../data/subj-{subj}_{attr}_rep-{rep}.csv")
    
    if rep > 0:
        fl = glob(f"../data/{subj}_{rep}-{attr}*mat")
        m = loadmat(fl[0])
        