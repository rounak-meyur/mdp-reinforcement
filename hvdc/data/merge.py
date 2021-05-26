# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 21:53:51 2021

Author: Rounak Meyur
Description: Merges the raw current and voltage time series data to construct 
the training dataset.
"""

import os
import pandas as pd



def merge_data(path):
    filelist = [f for f in os.listdir(path) if f.endswith(".csv")]
    # Get columns from individual csv files and merge them
    df_current = pd.DataFrame(columns=['time'])
    df_voltage = pd.DataFrame(columns=['time'])
    for f in filelist:
        df = pd.read_csv("raw/"+f)
        cname = '/'.join(f[4:-4].split('_'))
        if f.startswith("cur"):
            df_current[cname] = df[" Ia2"]
        elif f.startswith("vol"):
            df_voltage[cname] = df[" Vbus"]
    df_current['time'] = df["Domain"]
    df_voltage['time'] = df["Domain"]
    
    # Write to new csv file
    df_current.to_csv('current.csv',index=False)
    df_voltage.to_csv('voltage.csv',index=False)
    return


merge_data(path="raw/")
