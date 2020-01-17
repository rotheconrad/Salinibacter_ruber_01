#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -fad fastANI_directory -op output_prefix
## This script reads one vs all fastANI files from a directory and plots ANI vs Shared Fragments / Total Fragments.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Friday, June 21, 2019
## Date Updated :: N/A first version

def ANI_vs_shared_plot(df, org, op):
    ''' Coordinates reading fastANI directory files and plotting '''

    # Calculate Correlations
    pcorr = round(df['xs'].corr(df['ys'], method='pearson'), 4)
    scorr = round(df['xs'].corr(df['ys'], method='spearman'), 4)
    kcorr = round(df['xs'].corr(df['ys'], method='kendall'), 4)
    pcorrt = f'Pearson Correlation Coefficient: {pcorr}'
    scorrt = f"Spearman's Correlation: {scorr}"
    kcorrt = f"Kendall's Tau: {kcorr}"

    # Calculate mean and median
    ani_mean = round(np.mean(df['xs']), 2)
    ani_median = round(np.median(df['xs']), 2)

    #sns.set(style="white")
    plt.figure(figsize=(20, 10))
    ax = plt.gca()
    oName = ' '.join(org.split('_'))
    ax.text(0.5, 1.13, oName, fontstyle='italic', fontweight='heavy', color='#933b41', fontsize=60, horizontalalignment='center', transform=ax.transAxes)
    plt.title('ANI Value vs Shared Sequence Fragment Ratio', color='#933b41', fontsize=50, y=1.02)

    # Plot mean and median
    mn = plt.axvline(x=ani_mean, ymin=0, ymax=1, color='#525252', linewidth=4, linestyle=':', label='Mean')
    md = plt.axvline(x=ani_median, ymin=0, ymax=1, color='#737373', linewidth=4, linestyle='--', label='Median')

    # Write correlations to upper left corner
    ax.text(0.01, 0.99, f"{pcorrt}\n{scorrt}\n{kcorrt}", verticalalignment='top', fontsize=22, color='#737373', transform=ax.transAxes)

    # Build legend for mean and median
    plt.legend(handles=[mn, md], loc='lower right', fontsize=30, frameon=False, ncol=2)

    # Plot data points
    plt.scatter(df['xs'], df['ys'], marker='o', s=30, color="#933b41", alpha=0.15)

    # Setup labels, ticks, grid, etc
    plt.xlabel('ANI', fontsize=28, y=-0.02)
    plt.ylabel('Ratio (Shared Fragments / Total Fragments)', fontsize=28)
    ax.set_xticks(np.arange(94.9, 100.1, 0.1), minor=True)
    plt.xticks(np.arange(95.0, 100.1, 0.5), fontsize=18)
    ax.set_yticks(np.arange(0.45, 1.02, 0.02), minor=True)
    plt.yticks(np.arange(0.5, 1.01, 0.1), fontsize=18)
    plt.tick_params(axis='x', which='major', direction='in', color='k', width=6, length=12, bottom=True)
    ax.set_xlim(left=94.9, right=100.1)
    ax.set_ylim(bottom=0.46, top=1.02)
    ax.yaxis.grid(which="minor", color='#d9d9d9', linestyle='--', linewidth=1)
    ax.xaxis.grid(which="minor", color='#f0f0f0', linestyle='-', linewidth=1)
    ax.yaxis.grid(which="major", color='#d9d9d9', linestyle='--', linewidth=1.5)
    ax.xaxis.grid(which="major", color='#f0f0f0', linestyle='-', linewidth=2)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    plt.subplots_adjust(top=0.8)
    #plt.gcf().set_tight_layout(True)
    plt.savefig(f'{op}_plot.pdf')
    plt.close()

    with open(f'{op}_summary.tsv', 'w') as o:
        header=f'{org}\t{pcorr}\t{scorr}\t{kcorr}\n'
        xline='\t'.join(str(x) for x in df['xs'])
        yline='\t'.join(str(y) for y in df['ys'])
        o.write(header)
        o.write(f'{xline}\n')
        o.write(f'{yline}\n')

def gather_data(fad):
    ''' reads through the directory of fastANI file names, extracts the data, and returns a dataframe '''
    data_dict = {'xs': [], 'ys': []}
    f_list = os.listdir(fad)
    for file in f_list:
        with open(f'{fad}{file}', 'r') as f:
            for l in f:
                X = l.rstrip().split('\t')
                ani = float(X[2])
                shared = float(X[3])
                total = float(X[4])
                ratio = shared / total
                data_dict['xs'].append(ani)
                data_dict['ys'].append(ratio)

    df = pd.DataFrame(data_dict)
    df = df[df['xs'] != 100.0]
    return df

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(description='This script reads one vs all fastANI files from a directory and plots ANI vs Shared Fragments / Total Fragments.')
    parser.add_argument('-fad', '--fastani_directory', help='Please specify the directory with one vs all fastANI files!', required=True)
    parser.add_argument('-org', '--organism_name', help='What organism are these genomes from? Use underscore. Expl: Salinibacter_ruber', required=True)
    parser.add_argument('-op', '--output_prefix', help='What do you want to name the output file?', required=True)
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')
    df = gather_data(args['fastani_directory'])
    ANI_vs_shared_plot(df, args['organism_name'], args['output_prefix'])

if __name__ == "__main__":
    import argparse, os, matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
#    import seaborn as sns
    main()
