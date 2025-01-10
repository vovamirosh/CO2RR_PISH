import pandas as pd
import csv
import sqlite3
import os


def create_sq():
    dir_name = os.getcwd()

    df_to_db = pd.read_csv(f"{dir_name}\Datasets\data_with_0.csv", index_col=0)
    conn = sqlite3.connect(f'{dir_name}\Datasets\\final.db')

    df_to_db.to_sql('co2rr', conn, if_exists='replace', index=False)

    conn.close()