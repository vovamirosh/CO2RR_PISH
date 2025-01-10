import pandas as pd
import numpy as np
import os
from io import BytesIO

import requests
import pandas as pd


def preprocess(link):
    r = requests.get(link)
    data = r.content

    dir_name = os.getcwd()

    fe = pd.read_csv(BytesIO(data), index_col=0)
    fe.to_csv(f"{dir_name}\\Datasets\\raw_dataset.csv")
    
    fe_cl = fe.copy()

    #Delete columns as they are non-informative for our research
    fe_cl = fe_cl.drop(["Number", "Autor", "Publication", "Title", "Journal"], axis=1)

    fe_cl = fe_cl.dropna(subset=["Formula"])

    # Drop empty columns
    empt_col = fe_cl.columns[fe_cl.isna().all()]
    fe_cl = fe_cl.drop(empt_col, axis=1)

    #Delete columns because not of all articles have such parameter
    to_drop = ['Unnamed: 0', 'IF', 'Capacitance, (uF/cm2)', 'Current density, mA/cm2','ECSA, cm2', 'BET surface, m2/g', 'Loading, μg/cm2']
    fe_cl = fe_cl.drop([col_to_dr for col_to_dr in to_drop if col_to_dr in list(fe_cl.columns)], axis=1)

    #Drop row with filter from Google sheet
    fe_cl = fe_cl.drop(fe_cl.index[0])


    # 1. we can see that in support column there is duplicates with space. We need to rename them to a normal format
    # 2. we can see that in the df there are catalyst with complex structure such as: C0.902O0.074N0.024. It is important to drop them due to their class
    # 3. Also there is a duplication in a Base name. Some of them are one, but called in a different name

    fe_cl["Base"] = fe_cl["Base"].replace('gas diffusion layer', "gasdiffusion electrode")
    val_to_drop = fe_cl['Support'].isin(['C0.902O0.074N0.024', 'C0.912O0.054N0.017'])
    fe_cl = fe_cl[~val_to_drop]
    fe_cl = fe_cl.reset_index(drop=True)

    # Firstly I replaced nan values in the "Support" columns by string "No", after deleted additional spaces in the compounds' names
    fe_cl["Support"] = fe_cl["Support"].fillna(value="No")
    for i in range(0, len(fe_cl)):
        val = fe_cl.loc[i, "Cat. compound"]
        try:
            fe_cl.loc[i, "Cat. compound"] = val.replace(" ", "")
        except AttributeError:
            print(val)

    for i in range(0, len(fe_cl)):
        val = fe_cl.loc[i, "Support"]
        fe_cl.loc[i, "Support"] = val.replace(" ", "")

    # For some geometry there ia no length so imply nan like 0
    fe_cl["length min (nm)"] = fe_cl["length min (nm)"].fillna(0)
    fe_cl["length max (nm)"] = fe_cl["length max (nm)"].fillna(0)
    fe_cl["length aver (nm)"] = fe_cl["length aver (nm)"].fillna(0)
    fe_cl["Pore size, nm"] = fe_cl["Pore size, nm"].fillna(0)

    # For some catalyst test no information about time, but the most common is 1 hour
    fe_cl["Time, h"] = fe_cl["Time, h"].fillna(1)

    # Some catalysts concsist from one metall so there is no support compound
    fe_cl["Support"] = fe_cl["Support"].fillna("No")
    fe_cl["Cat./Support ratio, % (at)"] = fe_cl["Cat./Support ratio, % (at)"].fillna(100)

    # I am filling information about the average size of particles, there is no information about it
    nan_width = fe_cl[fe_cl["width aver (nm)"].isna()].index
    fe_cl.loc[nan_width, "width aver (nm)"] = (fe_cl.loc[nan_width,"width min (nm)"]+ fe_cl.loc[nan_width,"width max (nm)"])/2
    nan_len = fe_cl[fe_cl["length aver (nm)"].isna()].index
    fe_cl.loc[nan_len, "length aver (nm)"] = (fe_cl.loc[nan_len, "length min (nm)"] + fe_cl.loc[nan_len, "length max (nm)"])/2


    df_len = len(fe_cl)
    to_drop_prod = [prod for prod in fe_cl["Product"].unique() if len(fe_cl[fe_cl["Product"] == prod]) / df_len < 0.01]
    fe_cl= fe_cl[~fe_cl["Product"].isin(to_drop_prod)]

    print(to_drop_prod)

    fe_cl["Product"] = fe_cl["Product"].replace("CH24", "C2H4")

    print(fe_cl["Product"].unique())

    fe_cl.to_csv(f"{dir_name}\\Datasets\\prepr_.csv", index=False)

if __name__ == "__main__":
    dir_name = os.getcwd()
    preprocess('https://docs.google.com/spreadsheets/d/e/2PACX-1vQxIFJUE5JeauvxO11xI_LXezw_LNB_Rv4CNYoDJ0EIFKkLJiZ4ERt4CU4V5R2mEQvCp-n_A2OqKwx_/pub?gid=1053270247&single=true&output=csv')