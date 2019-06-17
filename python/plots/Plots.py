#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from IntaRNApvalue import IntaRNApvalue
matplotlib.use('Agg')  # run matplotlib headless


class Plots:
    @staticmethod
    def get_dist_function(query: str, target: str, n: int, shuffle_mode='b') -> None:
        """Gets a distribution function as bar histogram to a target/query sequence combination"""
        i = IntaRNApvalue(['-q', query, '-t', target, '-a', str(n), '-sm', shuffle_mode, '--threads', '0'])
        scores, non_interactions = i.get_scores()
        percent_non_int = round(non_interactions / (n + non_interactions) * 100, 1)
        annot_non_int = '{}% of all sequence pairs had no interaction'.format(percent_non_int)

        fig, ax = plt.subplots()
        ax.annotate(annot_non_int, (0.05, 0.01), rotation=90, size=10)
        n, bins, patches = ax.hist(scores, 100, density=True, facecolor='g', range=(min(scores), 0))
        plt.xlabel('MFE')
        plt.ylabel('MFE Frequency')
        plt.title('Distribution of IntaRNA scores from shuffled sequences')
        plt.grid(True)

        # Try to fit a gaussian distribution
        avg = np.mean(scores)  # average
        var = np.var(scores)  # variance
        std_dev = np.sqrt(var)  # standard deviation

        gauss_x = np.linspace(np.min(scores), np.max(scores), 100)  # even spaced numbers over interval
        gauss_y = 1.0 / np.sqrt(2 * np.pi * var) * np.exp(-0.5 * ((gauss_x - avg) ** 2 / var))
        # f(x, mu, std_dev) = 1/std_dev*sqrt(2*pi) * e^(-0.5 * ((x-mu)/std_dev)^2)
        plt.plot(gauss_x, gauss_y, 'k--', label='Gauss distribution')

        # TODO: Gumbel distribution

        plt.legend(loc='upper left')
        plt.savefig('dist_histogram_{}'.format(n))

    @staticmethod
    def get_pvalue_graph(query: str, target: str, max_exp: int = 4) -> None:
        """Plots the p-value for different shuffle types and # of scores, max_exp defines # scores with 10^max_exp
        For higher order of exponents this will take a very long time"""

        for n in range(1, max_exp + 1):
            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 'q', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'rx', label='only query shuffled' if n == 1 else '')

            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 't', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'bx', label='only target shuffled' if n == 1 else '')

            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 'b', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'gx', label='query & target shuffled' if n == 1 else '')

        plt.title('p-values for different shuffle types and shuffle amount')
        plt.xlabel('# of used scores')
        plt.ylabel('p-value')
        plt.legend(loc='upper left')
        plt.savefig('pvalue_graph_{}'.format(10**max_exp))


if __name__ == '__main__':
    q = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    t = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
        'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    # Plots.get_dist_function(q, t, 100, 'b')
    Plots.get_pvalue_graph(q, t, 4)
