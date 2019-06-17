#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

import matplotlib
matplotlib.use('Agg')  # run matplotlib headless
from matplotlib import pyplot as plt
from scipy.stats import genextreme as gev
import numpy as np

import sys
sys.path.insert(0, '..')  # add parent so we can import
from IntaRNApvalue import IntaRNApvalue


class Plots:
    @staticmethod
    def generate_histogram(query: str, target: str, shuffles: int, shuffle_mode='b', zoomed=False) -> None:
        """Gets a distribution function as bar histogram to a target/query sequence combination"""
        i = IntaRNApvalue(['-q', query, '-t', target, '-a', str(shuffles), '-sm', shuffle_mode, '--threads', '0'])
        scores, non_interactions = i.get_scores()
        scores = [-score for score in scores]  # invert for better fitting because of negative skew

        percent_non_int = round(non_interactions / (shuffles + non_interactions) * 100, 1)
        annot_non_int = '{}% of all sequence pairs had no interaction'.format(percent_non_int)

        fig, ax = plt.subplots(figsize=(11.69, 8.27))  # DIN A4 size
        ax.annotate(annot_non_int, (-0.35, 0.01), rotation=90, size=10)
        n, bins, patches = ax.hist(scores, 100, density=True, facecolor='g', range=(0, max(scores)))
        plt.xlabel('MFE (* -1)')
        plt.ylabel('Density')
        shuffle_mode_str = 'query only' if shuffle_mode == 'q' else 'target only' if shuffle_mode == 't' else 'both'
        plt.title('Distribution of IntaRNA scores from shuffled sequences\nShuffle mode: {}'.format(shuffle_mode_str))
        plt.grid(True)

        avg = np.mean(scores)  # average
        var = np.var(scores)  # variance
        std_dev = np.sqrt(var)  # standard deviation

        # GAUSSIAN DISTRIBUTION
        # f(x, mu, std_dev) = 1/std_dev*sqrt(2*pi) * e^(-0.5 * ((x-mu)/std_dev)^2)
        plt.plot(bins, 1.0 / np.sqrt(2 * np.pi * var) * np.exp(-0.5 * ((bins - avg) ** 2 / var)),
                 'k--', label='Gauss distribution')

        # GUMBEL DISTRIBUTION
        # mu is location parameter and beta is scale parameter
        # https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
        beta = std_dev * np.sqrt(6) / np.pi
        mu = avg - np.euler_gamma * beta
        plt.plot(bins, (1 / beta) * np.exp(-(bins - mu) / beta) * np.exp(-np.exp(-(bins - mu) / beta)),
                 'r--', label='Gumbel distribution')

        # GENERALIZED EXTREME VALUE DISTRIBUTION
        fit = gev.fit(scores)  # TODO: c, loc, scale parameters
        pdf = gev.pdf(bins, *fit)
        plt.plot(bins, pdf, 'b--', label='Generalized Extreme Value distribution')

        plt.legend(loc='upper left')
        plt.savefig('dist_histogram_sm={}_{}.png'.format(shuffle_mode, shuffles),
                    orientation='landscape', format='png', dpi=300)

    @staticmethod
    def generate_pvalue_graph(query: str, target: str, max_exp: int = 4) -> None:
        """Plots the p-value for different shuffle types and # of scores, max_exp defines # scores with 10^max_exp
        For higher order of exponents this will take a very long time"""

        plt.figure(figsize=(11.69, 8.27))  # DIN A4 size

        for n in range(1, max_exp + 1):
            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 'q', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'rx', label='only query shuffled' if n == 1 else '')

            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 't', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'bx', label='only target shuffled' if n == 1 else '')

            i = IntaRNApvalue(['-q', query, '-t', target, '--amount', str(10**n), '--shuffle', 'b', '--threads', '0'])
            plt.semilogx(10 ** n, i.calculate_pvalue_gauss(), 'gx', label='query & target shuffled' if n == 1 else '')

        plt.title('p-values for different shuffle modes and amount of used scores')
        plt.xlabel('# scores')
        plt.ylabel('p-value')
        plt.legend(loc='upper left')
        plt.savefig('pvalue_graph_{}.png'.format(10**max_exp),
                    orientation='landscape', format='png', dpi=300)


if __name__ == '__main__':
    q = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    t = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
        'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    Plots.generate_histogram(q, t, 1000, 'b')
    """
    Plots.generate_histogram(q, t, 1000, 'q')
    Plots.generate_histogram(q, t, 10000, 'q')
    Plots.generate_histogram(q, t, 100, 't')
    Plots.generate_histogram(q, t, 1000, 't')
    Plots.generate_histogram(q, t, 10000, 't')
    Plots.generate_histogram(q, t, 100, 'b')
    Plots.generate_histogram(q, t, 1000, 'b')
    Plots.generate_histogram(q, t, 10000, 'b')
    """
    # Plots.generate_pvalue_graph(q, t, 1)
