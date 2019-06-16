#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

from matplotlib import pyplot as plt
import numpy as np
from IntaRNApvalue import IntaRNApvalue


class Plots:
    def __init__(self):
        self.intarna = IntaRNApvalue()

    def get_dist_function(self) -> None:
        """Gets a distribution function as bar histogram to a target/query sequence combination"""
        scores = self.intarna.get_scores()

        fig, ax = plt.subplots()
        n, bins, patches = ax.hist(scores, 100, density=True, facecolor='g', range=(min(scores), 0))
        plt.xlabel('MFE')
        plt.ylabel('Wahrscheinlichkeit')
        plt.title('Wahrscheinlichkeit eines IntaRNA scores aus randomisierten Sequenzen')
        plt.grid(True)

        # Try to fit a gaussian distribution
        avg = np.mean(scores)  # average
        var = np.var(scores)  # variance
        std_dev = np.sqrt(var)  # standard deviation
        # gauss_x = np.linspace(np.min(scores), np.max(scores), 100)  # even spaced numbers over interval
        gauss_x = np.linspace(np.min(scores), 0, 100)  # even spaced numbers over interval
        gauss_y = 1.0 / np.sqrt(2 * np.pi * var) * np.exp(-0.5 * ((gauss_x - avg) ** 2 / var))
        # f(x, mu, std_dev) = 1/std_dev*sqrt(2*pi) * e^(-0.5 * ((x-mu)/std_dev)^2)
        plt.plot(gauss_x, gauss_y, 'k--')

        plt.show()

    def get_pvalue_graph(self, query: str, target: str, max_exp: int = 4) -> None:
        """Plots the p-value for different shuffle types and amounts, max_exp defines amount with 10^max_exp"""
        self.intarna.query = query
        self.intarna.target = target
        for n in range(1, max_exp):
            self.intarna.shuffle_query, self.intarna.shuffle_target = True, False
            plt.plot(10 ** n, self.intarna.calculate_pvalue_empirical(), 'rx')
            self.shuffle_query, self.shuffle_target = False, True
            plt.plot(10 ** n, self.intarna.calculate_pvalue_empirical(), 'bx')
            self.intarna.shuffle_query, self.intarna.shuffle_target = True, True
            plt.plot(10 ** n, self.intarna.calculate_pvalue_empirical(), 'gx')
        plt.title('p-values for different shuffle types and shuffle amount')
        plt.xlabel('# of shuffles')
        plt.ylabel('p-value')
        plt.show()


if __name__ == '__main__':
    q = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    t = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
        'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    p = Plots()
    p.get_dist_function()
    # p.get_pvalue_graph(q, t, 4)
