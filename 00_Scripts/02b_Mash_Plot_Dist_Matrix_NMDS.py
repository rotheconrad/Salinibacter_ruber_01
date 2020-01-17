#!/usr/bin/env python

'''NMDS Plot for Mash Distance Matrix

Writes plot in .png format

# The stress value is calculated, normalized and printed on the plot.
# The stress value can be interpretted as follows:
# According to Kruskal (1964, p. 3): value 0 indicates "perfect" fit, 0.025 excellent,
# 0.05 good, 0.1 fair, and 0.2 poor.
# For more information cf. Kruskal (1964, p. 8-9) and Borg (2005, p. 41-43).

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: November 22nd, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.manifold import MDS

def NMDS_Plot(df, outfile):

        m = '.'

        c = '#ACDE64'

        seed = np.random.RandomState(seed=3)
        nmds = MDS(n_components=2, metric=False, max_iter=10000, eps=1e-9,
                dissimilarity="precomputed", random_state=seed, n_jobs=1, n_init=10000)

        npos = nmds.fit_transform(df)

        distances = []

        for i in range(len(npos)-1):
                j = i+1
                distances.append(np.sqrt( (npos[j,0]-npos[i,0])**2 + (npos[j,1]-npos[i,1])**2) )

        stress = np.sqrt(nmds.stress_ / sum(distances)**2)
        stext = 'Stress = %f' % (stress)

        for i in range(len(npos)):
                plt.scatter(npos[i, 0], npos[i, 1], c=c, marker=m, label=df.index[i])

        plt.subplots_adjust(right=0.7)
        ax = plt.gca()
        plt.axis('equal')
        plt.title('NMDS plot of Mash Distance')
        #plt.legend(frameon=False, bbox_to_anchor=(1.04,1.05), loc="upper left")
        plt.text(
            1, .010, stext, fontsize=10, color='#737373',
            horizontalalignment='right', transform=ax.transAxes
            )

        plt.savefig(outfile)
        plt.close()         

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the Mash Distance Matrix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='Please specify the output file name! (name.png)',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    print('\nRunning Script...\n')

    df = pd.read_csv(args['input_file'], sep='\t', index_col=0)

    _ = NMDS_Plot(df, args['out_file'])

    print('\nLooks like the plot was created successfully!\n\n')

if __name__ == "__main__":
    main()
