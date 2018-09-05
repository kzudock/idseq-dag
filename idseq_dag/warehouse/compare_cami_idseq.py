#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import idseq_dag.util.command as command
import sample_lists
import data_wrangling as dw

def main():
    warehouse_dir = "/mnt/idseq/warehouse"
    figure_dir = f"{warehouse_dir}/figures"
    command.execute(f"mkdir -p {figure_dir}")
    idseq_predictions = f"{warehouse_dir}/idseq_predictions.csv"
    cami_labels = f"{warehouse_dir}/cami_labels.csv"
    idseq_df = pd.read_csv(idseq_predictions, index_col=[0])
    cami_df = pd.read_csv(cami_labels, index_col=[0])
    for sample_name in sample_lists.CAMI_Airways_sample_names:
        clean_sample_name = sample_lists.clean_name(sample_name)
        idseq_genus = dw.subset(idseq_df, tax_level=2, clean_sample_name)
        cami_genus = dw.subset(cami_df, tax_level=2, clean_sample_name)
        df = pd.merge(idseq_genus, cami_genus, how='outer', on='taxid', suffixes=('_idseq', '_cami'))
        df = df.fillna(0)
        rms = np.sqrt(mean_squared_error(df['count_cami'], df['count_idseq']))
        plt.ticklabel_format(style='sci', scilimits=(0,0))
        plt.scatter('count_idseq', 'count_cami', data=df, c='b')
        plt.xlabel('idseq count ((NT.r + NR.r)/2)')
        plt.ylabel('cami count')
        plt.title(f"{sample_name}\n(RMS error = {round(rms, 2)})")
        plt.savefig(f"{figure_dir}/{clean_sample_name}.png", format="png")
    command.execute(f"aws s3 cp {figure_dir} s3://idseq-samples-test/charles/ --recursive")

if __name__ == "__main__":
    main()
