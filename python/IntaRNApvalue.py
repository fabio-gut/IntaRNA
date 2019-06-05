#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

import os
import sys
import argparse
import multiprocessing as mp
from typing import Tuple, List
import time
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import quad as integ
from DinuclShuffle import dinucl_shuffle as din_s


class IntaRNApvalue:
    def __init__(self):
        self.bin = self.find_binary()

    @staticmethod
    def find_binary() -> str:
        """Tries to find the IntaRNA executable and returns path to binary or exits with error code 1 if not found"""
        to_search = os.environ['PATH'].split(':')  # get path env variable
        to_search.append(os.path.abspath(os.path.join(os.curdir, '../..')))  # add 2 dirs above this script
        bin_name = 'IntaRNA.exe' if os.name == 'nt' else 'IntaRNA'
        for path in to_search:  # search in each directory including subdirectories
            for dir_path, dir_name, file_names in os.walk(path):
                if bin_name in file_names:
                    return os.path.join(dir_path, bin_name)
        else:
            print('Error: Cannot find IntaRNA binary executable, please add it to your PATH')
            sys.exit(1)

    def process_cmd_args(self) -> None:
        """Processes all commandline args"""
        parser = argparse.ArgumentParser(description='Calculates p-values to IntaRNA scores')
        parser.add_argument('--query', type=str, required=True, help='Query sequence')
        parser.add_argument('--target', type=str, required=True, help='Target sequence')
        parser.add_argument('--amount', type=int, required=True, help='How often he sequences are shuffled')
        parser.add_argument('--shuffle-query', type=bool, default=True, dest='sq',
                            help='Amount of shuffles to the query sequence')
        parser.add_argument('--shuffle-target', type=int, default=True, dest='st',
                            help='Amount of shuffles to the target sequence')
        args = parser.parse_args()
        self.calculate_pvalue_empirical(args.query, args.target, args.amount, args.sq, args.st)

    def get_score(self, seq_tuple: Tuple[str, str]) -> float:
        """Gets the IntaRNA score of a single query/target tuple, query is first element, target second"""
        res = os.popen('{} -q {} -t {}'.format(self.bin, seq_tuple[0], seq_tuple[1])).read()
        if 'no favorable interaction for target and query' in res:
            return 0.1  # TODO: What do with query/target combinations that don't interact?
        return float(res.split('interaction energy = ')[1].split(' kcal/mol\n')[0])

    def get_scores(self, seq_list: List[Tuple[str, str]]) -> List[float]:
        """Gets the IntaRNA score to a given list of query/target tuples, returns a sorted list of scores"""
        pool = mp.Pool(processes=mp.cpu_count())  # we will do this in parallel and not sequential
        scores = pool.map(self.get_score, seq_list)
        scores.sort()  # Timsort is O(n log n) on average and worst, O(n) on best case
        return scores

    @staticmethod
    def shuffle_sequence(query: str, target: str, n: int, shuffle_q: bool, shuffle_t: bool) -> List[Tuple[str, str]]:
        """Shuffles a query/target pair n times and returns an array of sequence tuples, by default both are shuffled.
            The returned array always has length n+1, duplicate entries are possible"""
        shuffles = [(query, target)]  # TODO: add original seq???
        for _ in range(n):
            shuffles.append((din_s(query) if shuffle_q else query, din_s(target) if shuffle_t else target))
        return shuffles

    def calculate_pvalue_empirical(self, query: str, target: str, n: int, shuffle_q: bool, shuffle_t: bool) -> float:
        """Calculates a p-value to a target/query combination empirical with a given amount of shuffle iterations"""
        original_score = self.get_score((query, target))
        shuffles = self.shuffle_sequence(query, target, n, shuffle_q, shuffle_t)
        scores = self.get_scores(shuffles)
        return [score <= original_score for score in scores].count(True) / len(scores)

    def calculate_pvalue(self, query: str, target: str, n: int, shuffle_q: bool, shuffle_t: bool) -> any:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations"""
        original_score = self.get_score((query, target))
        shuffles = self.shuffle_sequence(query, target, n, shuffle_q, shuffle_t)
        scores = self.get_scores(shuffles)

        # Try to fit a gaussian distribution
        avg = np.mean(scores)  # average
        var = np.var(scores)  # variance
        std_dev = np.sqrt(var)  # standard deviation

        def gauss(x):
            return 1.0 / np.sqrt(2 * np.pi * var) * np.exp(-0.5 * ((x - avg) ** 2 / var))
        return integ(gauss, -np.inf, original_score)

    def get_dist_function(self, query: str, target: str, n: int, shuffle_q: bool, shuffle_t: bool) -> None:
        """Gets a distribution function as bar histogram to a target/query sequence combination"""
        shuffles = self.shuffle_sequence(query, target, n, shuffle_q, shuffle_t)
        scores = self.get_scores(shuffles)

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
        """Plots the p-value for different shuffle types and amounts, max_exp defines amount with 10^max_exp
        TODO: something is weird, doesnt plot all values?!?!"""
        for n in range(1, max_exp):
            plt.plot(10 ** n, self.calculate_pvalue_empirical(query, target, 10 ** n, True, False), 'rx')
            plt.plot(10 ** n, self.calculate_pvalue_empirical(query, target, 10 ** n, False, True), 'bx')
            plt.plot(10 ** n, self.calculate_pvalue_empirical(query, target, 10 ** n, True, True), 'gx')
        plt.title('p-values for different shuffle types and shuffle amount')
        plt.xlabel('# of shuffles')
        plt.ylabel('p-value')
        plt.show()


if __name__ == '__main__':
    q = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    t = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
        'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    i = IntaRNApvalue()
    # i.process_cmd_args()
    i.get_dist_function(q, t, 10000, True, True)
    # i.get_pvalue_graph(q, t, 4)
    # print('Integriert: {}'.format(i.calculate_pvalue(q, t, 1000, True, True)))
    # print('Empirisch:  {}'.format(i.calculate_pvalue_empirical(q, t, 1000, True, True)))
