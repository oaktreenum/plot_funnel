# helical_bundle_predict and fast_relax_near_native out log file to funnel plot script, 
# Allon Goldberg, Research Assistant, Flatiron Institute, 2/2025

# Usage: python plot_FR+HBP_funnel.py <any_string>
# [SEE TREE BELOW] Run in dir of `bundleout.log`, with `score.sc` in FRfill/

# RUN_DIRECTORY/
# ├── bundleout.log
# └── FRfill/
#     └── score.sc


import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
from scipy.special import logsumexp

current_dir = os.path.basename(os.getcwd())

# Helper function: extract helical bundle predict output from the out file named bundleout.log
def hbp_extract(log_files):
    # Gather input files into vector according to given pattern
    log_names = glob.glob(log_files)

    for file_name in log_names:
        file_name_sep = file_name.split('.')
        # Read file
        with open(file_name, 'r') as log:
            lines = log.readlines()

        p_near = 0
        Qp_near = 0
        df_lines = []
        mode = 0
        for line in lines:
            # Find where header starts and change mode to read data
            if 'MPI_worker_node' in line:
                mode = 1
                df_lines.append(line)
                continue
            # If end of data, change mode
            if 'End' in line:
                mode = 0
            # Save data in vector to turn into data frame later
            if mode == 1:
                df_lines.append(line)
                continue
            # Find PNear value
            if 'PNear:' in line:
                split = line.split()
                p_near = float(split[1])
                continue
            if 'PNear_additional_1_weights/qm_GAMESS.wts:' in line:
                Qsplit = line.split()
                Qp_near = float(Qsplit[1])
                break
            
        fake_file = '\n'.join(df_lines)
        data_io = StringIO(fake_file)
        df_raw = pd.read_csv(data_io, sep='\t{1,2}', engine='python', index_col=False).fillna(0)

        # Minimum Energies for Axes Scaling
        Energy_min = df_raw['Energy'].min()
        QEnergy_min = df_raw['Additional_score_1_weights/qm_GAMESS.wts'].min()
        QEnergy_2ndmin = df_raw[df_raw['Additional_score_1_weights/qm_GAMESS.wts'] > QEnergy_min]['Additional_score_1_weights/qm_GAMESS.wts'].min()

        #df = df_raw[(df_raw['RMSD'] <= RMSD_98th) & (df_raw['Energy'] <= Energy_998th) & (df_raw['RMSD_to_best_for_score_1_weights/qm_GAMESS.wts'] <= QRMSD_98th) & (df_raw['Additional_score_1_weights/qm_GAMESS.wts'] <= QEnergy_90th)]
        df = df_raw


    return df, p_near, Qp_near, Energy_min, QEnergy_min, QEnergy_2ndmin
        # fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        # classical_scatter_plot = sns.scatterplot(data=df[df['Energy'] <= Energy_min + 340], x='RMSD_to_best', y='Energy', hue="Hbonds", palette="crest", ax=axes[0])
        # axes[0].set_title('Classical Energy Calculation')
        # axes[0].annotate(f'PNear: {p_near}', xy=(0.6, 0.05), xycoords='axes fraction', fontsize=13, color='black')
        # axes[0].set_xlabel('RMSD')
        # axes[0].set_ylabel('Energy')
        # # QM_scatter_plot = sns.scatterplot(data=df[df['Additional_score_1_weights/qm_GAMESS.wts'] <= QEnergy_min + 340], x="RMSD_to_best_for_score_1_weights/qm_GAMESS.wts", y="Additional_score_1_weights/qm_GAMESS.wts", hue="Hbonds", palette="ch:s=-.2,r=.6", ax=axes[1])
        # QM_scatter_plot = sns.scatterplot(data=df, x="RMSD_to_best_for_score_1_weights/qm_GAMESS.wts", y="Additional_score_1_weights/qm_GAMESS.wts", hue="Hbonds", palette="ch:s=-.2,r=.6", ax=axes[1])       
        # axes[1].set_title('Quantum-Chemical Energy Calculation')
        # axes[1].annotate(f'PNear: {Qp_near}', xy=(0.6, 0.05), xycoords='axes fraction', fontsize=13, color='black')
        # axes[1].set_xlabel('RMSD')
        # axes[1].set_ylabel('Energy')

        # # Adjust layout
        # #plt.tight_layout()
        # # Save the figure as a PNG file
        # plt.savefig(f'plots_{file_name_sep[0]}.png', format='png', dpi=300)

# Helper function: PNear
def Pnear(df1,df2,score_column):

    df1.rename(columns={'rmsd': 'RMSD', 'Additional_score_1_weights/qm_GAMESS.wts': 'QMEnergy'}, inplace=True)
    df2.rename(columns={'rmsd': 'RMSD', 'score': 'Energy', 'total_energy': 'QMEnergy'}, inplace=True)
    # Combine DataFrames for unified processing
    combined_df = pd.concat([df1, df2], ignore_index=True)
    
    # Constants
    lambda_val = 2.5 
    kBT = 0.62   

    # Subtract the minimum score for normalization
    combined_df['adjusted_score'] = combined_df[score_column] - combined_df[score_column].min()

    # Compute the numerator terms: -rmsd^2 / lambda^2 - score / kBT
    combined_df['numerator_exp'] = (
        -combined_df['RMSD']**2 / lambda_val**2
        - combined_df['adjusted_score'] / kBT
    )

    # Compute the denominator terms: -score / kBT
    combined_df['denominator_exp'] = -combined_df['adjusted_score'] / kBT

    # Use log-sum-exp for numerical stability
    numerator_lse = logsumexp(combined_df['numerator_exp'])
    denominator_lse = logsumexp(combined_df['denominator_exp'])

    # Compute Pnear
    log_pnear = numerator_lse - denominator_lse
    pnear = np.exp(log_pnear)
    return pnear

# Make plots and save
def main():
    df, p_near, Qp_near, Energy_min, QEnergy_min, QEnergy_2ndmin = hbp_extract("bundleout.log")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'{current_dir}')
    classical_scatter_plot = sns.scatterplot(data=df[df['Energy'] <= Energy_min + 134], x='RMSD', y='Energy', hue="Hbonds", palette="crest", ax=axes[0])
    axes[0].set_title('Classical Energy Calculation')
    axes[0].annotate(f'PNear: {p_near}', xy=(0.6, 0.05), xycoords='axes fraction', fontsize=13, color='black')
    axes[0].set_xlabel('RMSD')
    axes[0].set_ylabel('Energy')
    # QM_scatter_plot = sns.scatterplot(data=df[df['Additional_score_1_weights/qm_GAMESS.wts'] <= QEnergy_min + 340], x="RMSD_to_best_for_score_1_weights/qm_GAMESS.wts", y="Additional_score_1_weights/qm_GAMESS.wts", hue="Hbonds", palette="ch:s=-.2,r=.6", ax=axes[1])
    QM_scatter_plot = sns.scatterplot(data=df[(df['Additional_score_1_weights/qm_GAMESS.wts'] <= QEnergy_2ndmin + 134)], x="RMSD", y="Additional_score_1_weights/qm_GAMESS.wts", hue="Hbonds", palette="ch:s=-.2,r=.6", ax=axes[1])       
    axes[1].set_title('Quantum-Chemical Energy Calculation')
    axes[1].annotate(f'PNear: {Qp_near}', xy=(0.6, 0.05), xycoords='axes fraction', fontsize=13, color='black')
    axes[1].set_xlabel('RMSD')
    axes[1].set_ylabel('Energy')


    with open("FRfill/score.sc", 'r') as log:
        lines = log.readlines()

        df_lines = []
        mode = 0
        for line in lines:
            # Find where header starts and change mode to read data
            if 'SCORE:' in line:
                mode = 1
                df_lines.append(line)
                continue
            # If end of data, change mode
            if 'End' in line:
                mode = 0
            # Save data in vector to turn into data frame later
            if mode == 1:
                df_lines.append(line)
                continue
        
        FRfake_file = '\n'.join(df_lines)
        FRdata_io = StringIO(FRfake_file)
        FRdf_raw = pd.read_csv(FRdata_io, sep='\s+', index_col=False).fillna(0)

    # Minimum Energies for Axes Scaling
    Energy_min1 = FRdf_raw['score'].min()
    RMSD_min1 = FRdf_raw['rmsd'].min()
    QEnergy_min1 = FRdf_raw['total_energy'].min()
    pnear1 = Pnear(df,FRdf_raw,'Energy')
    qpnear1 = Pnear(df,FRdf_raw,'QMEnergy')
    axes[0].annotate(f'PNear: {pnear1:.8f}', xy=(0.6, 0.1), xycoords='axes fraction', fontsize=13, color='orange')
    axes[1].annotate(f'PNear: {qpnear1:.8f}', xy=(0.6, 0.1), xycoords='axes fraction', fontsize=13, color='orange')
    classical_scatter_plot1 = sns.scatterplot(data=FRdf_raw[FRdf_raw['Energy'] <= Energy_min + 134], x='RMSD', y='Energy', ax=axes[0])
    QM_scatter_plot1 = sns.scatterplot(data=FRdf_raw[FRdf_raw['QMEnergy'] <= QEnergy_2ndmin + 134], x='RMSD', y='QMEnergy', ax=axes[1])


    
    fig.savefig(f'filled_{current_dir}.png', dpi=300)
    fig.savefig(f'filled_{current_dir}.svg')

    return pnear1, qpnear1   




if __name__ == "__main__":
    import sys
    if sys.argv[1] == 'help':
        print("""Usage: python plot_FR+HBP_funnel.py go
                [SEE TREE BELOW] Run in dir of `bundleout.log`, with `score.sc` in FRfill/
                RUN_DIRECTORY(NAME_OF_DESIGN)/
                ├── bundleout.log
                └── FRfill/
                      └── score.sc""")
    else:
        logs = sys.argv[1]
        print(f'\n♡♡♡ Thank you for choosing me as your plotter. I work hard to plot your figures ♡♡♡ \n\nPLOTTING {current_dir}...\n')
        pnear1, qpnear1 = main()
        print(f'''Pnear for {current_dir}: 
                        PNear: {pnear1:.4f} 
                        QM-PNear {qpnear1:.4f}\n''')
        print(f'The plots for {current_dir} were generated and saved as png(300dpi)/svg files. ☺ ENJOY ☺\n')
        
