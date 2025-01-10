import pandas as pd
import numpy as np
import pymatgen.core as pmg

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pubchempy as pcp

def descr_adding():
    df_main = pd.read_csv("Datasets\prepr_.csv")

    compounds = set()

    compounds = compounds.union(set(df_main["Cat. compound"].unique()))
    compounds = compounds.union(set(df_main["Support"].unique()))
    compounds.remove("No")
    compounds = list(compounds)

    ion_df = pymat_desc(df_main, compounds)
    rdkit_df = rdkit_desc(df_main, compounds)
    df_with_av = save_df(df_main, rdkit_df, ion_df)
    save_df_with_null(df_with_av)
    
def pymat_desc(df_main, products):
    # Now we will use library pymatgen to get additional descriptors for catalysts
    ion_energy = dict()

    ion_energy["Formulas"] = products

    # Now manually I add 2 rows, that corresponds to:
    # 1. Main metall in the compound
    # 2. Oxidation state of this metall
    ion_energy["Free form"] = ['Cu', 'Pd', 'N', 'Cu', 'Cu', 'Cu', 'Au', 'Cu', 'Ag']
    ion_energy["Oxidation state"] = [1, 0, 0, 2, 0, 2, 0, 1, 0]

    # Now I am adding information about the energy of ionization of corresponsing element
    ion_energy["Energy of +1 ion"] = [pmg.Element(ion_energy["Free form"][i]).ionization_energies[ion_energy["Oxidation state"][i]] for i in range(0, len(ion_energy["Free form"]))]

    # Calculating an ion radius, if the element isn't in oxidation state 0, otherwise calcualte atomic radius
    ion_energy["Radius"] = []
    for i in range(0, len(ion_energy["Energy of +1 ion"])):
        ox_state = ion_energy["Oxidation state"][i]
        el = ion_energy["Free form"][i]
        if ion_energy["Oxidation state"][i] == 0:
            ion_energy["Radius"].append(pmg.Element(el).atomic_radius)
        else:
            ion_energy["Radius"].append(pmg.Element(el).ionic_radii[ox_state])   

    dict_ions = dict()
    for i in range(0, len(ion_energy["Formulas"])):
        dict_ions[ion_energy["Formulas"][i]] = {"Energy of +1 ion": ion_energy["Energy of +1 ion"][i],
                                                "Radius": ion_energy["Radius"][i]}
    dict_ions["No"] = {"Energy of +1 ion": 0,
                        "Radius": 0}

    # Now we formed a dataset, which rows are taken from the original dataset
    ion_df = pd.DataFrame()
    for catal in ["Cat. compound", "Support"]:
        for_df = {catal + "_Energy of +1 ion": [],
                catal + "_Metal_radius": []}
        for i in range(0, len(df_main)):
            comp = df_main.at[i, catal]
            for_df[catal + "_Energy of +1 ion"].append(dict_ions[comp]["Energy of +1 ion"])
            for_df[catal + "_Metal_radius"].append(dict_ions[comp]["Radius"])

        new = pd.DataFrame(data = for_df)
        print(new.columns)
        if len(ion_df) == 0:
            ion_df = new.copy()
        else:
            ion_df = pd.concat([ion_df, new], axis=1, join="inner")
    return ion_df

def corr_col(df, val):
    corr_matrix = df.corr(method="spearman").abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > val)]
    return to_drop

def rdkit_desc(df_main, products):
    CIDa = []
    compounds = list(products)
    cid_dict = dict()
    for cmp in compounds:
        for compound in pcp.get_compounds(cmp, 'name'):
            CIDa.append(compound.cid)
  
            cid_dict[cmp] = compound.cid
    cid_dict["Cu3N"] = 56841037

    comp_smiles = dict()
    for comp, cid in cid_dict.items():
        c = pcp.Compound.from_cid(cid)
        comp_smiles[comp] = c.canonical_smiles

    descriptor_names = list(rdMolDescriptors.Properties.GetAvailableProperties())
    get_descriptors = rdMolDescriptors.Properties(descriptor_names)
    num_descriptors = len(descriptor_names)

    descriptors_set = np.empty((0, num_descriptors), float)

    for _, value in comp_smiles.items():
        molecule = Chem.MolFromSmiles(value)
        descriptors = np.array(get_descriptors.ComputeProperties(molecule)).reshape((1,num_descriptors))
        descriptors_set = np.append(descriptors_set, descriptors, axis=0)
    
    df_comp = pd.Series(data=comp_smiles.keys())

    df_rdkit = pd.DataFrame(descriptors_set, columns= descriptor_names)
    combinedd = pd.concat([df_comp, df_rdkit], axis =1, join='inner')
    combinedd = combinedd.loc[:, combinedd.nunique() > 1]

    # It is so many descriptors so we need to drop high-correlated 
    num_col = combinedd.select_dtypes(include=["float64", "int64"]).columns.to_list()
    comb_low_cor = combinedd.drop(corr_col(combinedd[num_col], 0.9), axis=1).copy()
    num_col = comb_low_cor.select_dtypes(include=["float64", "int64"]).columns.to_list()

    # Now we need to fill our main dataset by the descriptors for 2 columns: "Support", "Cat. compound". 
    # To fastem algorithm we will form a dictionary, which will provide descriptors by formula

    combinedd = comb_low_cor.copy()
    rdkit_dict = dict()
    for i in range(0, len(combinedd)):
        rdkit_dict[combinedd.loc[i, 0]] = dict()
        for col in combinedd.columns[1:]:
            rdkit_dict[combinedd.loc[i, 0]][col] = combinedd.loc[i, col]
    rdkit_dict["No"] = {col: 0 for col in combinedd.columns[1:]}

    for catal in ["Cat. compound", "Support"]:
        for_df = {catal+"_"+col: [] for col in combinedd.columns[1:]}

    rdkit_df = pd.DataFrame()
    for catal in ["Cat. compound", "Support"]:
        for_df = {catal+"_"+col: [] for col in combinedd.columns[1:]}

        for i in range(0, len(df_main)):
            comp = df_main.at[i, catal]
            for col in combinedd.columns[1:]:
                for_df[catal+"_"+col].append(rdkit_dict[comp][col])

        new = pd.DataFrame(data = for_df)
  
        if len(rdkit_df) == 0:
            rdkit_df = new.copy()
        else:
            rdkit_df = pd.concat([rdkit_df, new], axis=1, join="inner")
        
    return rdkit_df

def save_df(df_main, rdkit_df, ion_df):
    df_for_save = pd.concat([df_main, rdkit_df, ion_df], axis=1, join="inner")

    # Lets merge feature columns of Cat. compound and support. That will help us to reduce the deminsionality of data
    # Some columns will be merged to average arifmetic, another part as harmonic average

    to_arifm_av = ["exactmw",	"lipinskiHBA",	"lipinskiHBD",	"NumHBA",	"NumHeavyAtoms"]
    to_harm_av = ["labuteASA",	"CrippenClogP",	"kappa1",	"kappa3", 'Energy of +1 ion', 'Metal_radius']

    ratio = df_for_save["Cat./Support ratio, % (at)"]/100
    df_with_av = df_for_save.copy()

    for arif in to_arifm_av:
        df_with_av["Av_ar_" + arif] = df_with_av["Cat. compound_"+arif] * ratio + df_with_av["Support_"+arif] * (1 - ratio)
        df_with_av = df_with_av.drop(["Cat. compound_"+arif, "Support_"+arif], axis=1)

    for harm in to_harm_av:
        df_with_av["Av_harm_" + harm] = 1 / (ratio / df_with_av["Cat. compound_" + harm] + (1 - ratio) / df_with_av["Support_"+harm])
        df_with_av["Av_harm_" + harm] = df_with_av["Av_harm_" + harm].fillna(df_with_av["Cat. compound_" + harm])
        df_with_av = df_with_av.drop(["Cat. compound_" + harm, "Support_" + harm], axis=1)

    df_with_av.to_csv("Datasets\\with_av_features.csv", index=False)
    return df_with_av

def save_df_with_null(df_with_av):
    unic = df_with_av.drop(["Product", "FE, %"], axis=1).copy()
    unic = unic.drop_duplicates()

    products = list(df_with_av["Product"].unique())
    exp_count = len(unic)
    prod_col = [[prod] * exp_count for prod in products]
    prod_col = list(np.array(prod_col).reshape(1, -1)[0])

    unic_mult = pd.concat([unic] * len(products)).reset_index(drop=True).copy()
    unic_mult["Product"] = prod_col

    merged_df = pd.merge(unic_mult, df_with_av, how='left')
    merged_df = merged_df.fillna(0)

    merged_df.to_csv("Datasets\\data_with_0.csv")

if __name__ == "__main__":
    descr_adding()