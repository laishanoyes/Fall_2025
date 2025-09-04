#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 14:08:22 2025

@author: kgjohnst
"""

import pandas as pd
import numpy as np
from collections import Counter as cnt


metadata=pd.read_csv('~/Downloads/sequences.csv')

usa_indexes=[k for k in range(metadata.shape[0]) if 'USA' in metadata['Geo_Location'][k] or 'usa' in metadata['Geo_Location'][k]]
usa_indexes=np.array(usa_indexes)

metadata=metadata.iloc[usa_indexes]
pangos=cnt(metadata['Pangolin'])

threshold=276

keys=list(pangos.keys())
variants=[keys[k] for k in range(len(keys)) if pangos[keys[k]]>=threshold]

accession_numbers=[]
for var in variants:
    temp=metadata.loc[metadata['Pangolin']==var]
    acc=metadata['Accession'].values
    acc_number=np.random.choice(acc,1)
    accession_numbers.append(acc_number[0])
