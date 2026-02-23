import numpy as np
import os
from medflow import train_med, sim_med
import pandas as pd

df2 = pd.read_csv("../output/compas.csv")
df2["c_charge_degree"] = (df2["c_charge_degree"] == "F").astype(int)

df2.to_csv("../output/compas_cleaned.csv", index=False)

base_path = './'  # Define the base path for file operations.
folder = '../output/'  # Define the folder where files will be stored.
path = os.path.join(base_path, folder, '')  # Combines the base path and folder into a complete path.
dataset_name = 'compas_cleaned'  # Define the name of the dataset.

if not (os.path.isdir(path)):  # checks if a directory with the name 'path' exists.
    os.makedirs(path)  # if not, creates a new directory with this name. This is where the logs and model weights will be saved.
    
## MODEL TRAINING
train_med(path=path, dataset_name=dataset_name, treatment='race', confounder=[], mediator=["age", "priors_count","c_charge_degree"], 
          outcome='is_recid', test_size=0.2, cat_var=["c_charge_degree"], sens_corr=None, seed_split=1,
          model_name=path + 'seed_1', trn_batch_size=64, val_batch_size=512, learning_rate=1e-4, seed=1,
          nb_epoch=3000, nb_estop=20, val_freq=1, emb_net=[164, 48, 32], int_net=[32, 24, 16])

## EFFECT ESTIMATION
# Path-specific effects
sim_med(path=path, dataset_name=dataset_name, model_name=path + 'seed_1', n_mce_samples=10000, seed=1, inv_datafile_name='1_path_100k_v2',
        cat_list=[0, 1], moderator=None)
