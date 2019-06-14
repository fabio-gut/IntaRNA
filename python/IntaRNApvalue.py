#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

from subprocess import PIPE, Popen, run
import random
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
        self.query = ''
        self.target = ''
        self.n = 0
        self.shuffle_query = False
        self.shuffle_target = False
        self.threads = ''

    @staticmethod
    def find_binary() -> str:
        """Tries to find the IntaRNA executable and returns path to binary or exits with error code 1 if not found

        >>> IntaRNApvalue.find_binary()
        'IntaRNA'
        """
        if not run('IntaRNA --version', shell=True, stdout=PIPE, stderr=PIPE).returncode:
            return 'IntaRNA'

        # if binary not found in path, search in parent of this script recursively
        bin_name = 'IntaRNA.exe' if os.name == 'nt' else 'IntaRNA'
        for dir_path, dir_name, file_names in os.walk(os.path.abspath(os.path.join(os.curdir, '..'))):
            if bin_name in file_names:
                return os.path.join(dir_path, bin_name)

        print('Error: Cannot find IntaRNA binary executable, please add it to your PATH')
        sys.exit(1)

    def process_cmd_args(self, test_args=None) -> None:
        """Processes all commandline args

        >>> i = IntaRNApvalue()
        >>> i.process_cmd_args(['-q', 'AGCGU', '-t', 'AAAGGCC', '--amount', '10', '--shuffle', 'b', '--threads', '3'])
        >>> i.query
        'AGCGU'
        >>> i.target
        'AAAGGCC'
        >>> i.n
        10
        >>> i.shuffle_query
        True
        >>> i.shuffle_target
        True
        >>> i.threads
        '3'
        """
        parser = argparse.ArgumentParser(description='Calculates p-values to IntaRNA scores')
        parser.add_argument('-q', '--query', dest='query', type=str, help='Query sequence', required=True)
        parser.add_argument('-t', '--target', dest='target', type=str, help='Target sequence', required=True)
        parser.add_argument('-a', '--amount', dest='amount', type=int, required=True,
                            help='How many randomly generated scores are used to calculate the p-value')
        parser.add_argument('-s', '--shuffle', dest='shuffle', required=True, choices=['q', 't', 'b'],
                            help='Which sequences are going to be shuffled: both, query only or target only')
        parser.add_argument('--threads', type=str, default='0', help='Sets the amount of threads used for IntaRNA')
        parser.add_argument('--seed', type=str, default=None,
                            help='Random seed to make sequence generation deterministic')

        args = parser.parse_args(test_args)
        shuffle_query = True if args.shuffle in ['b', 'q'] else False
        shuffle_target = True if args.shuffle in ['b', 't'] else False
        self.query, self.target, self.n = args.query, args.target, args.amount
        self.shuffle_query, self.shuffle_target, self.threads = shuffle_query, shuffle_target, args.threads
        random.seed(a=args.seed)

    @staticmethod
    def shuffle_sequence(seq: str, n: int) -> List[str]:
        """Shuffles a sequence n times and returns a list of sequences, duplicate entries are possible

        >>> random.seed('IntaRNA')
        >>> IntaRNApvalue.shuffle_sequence('AGGAUGGGGGA', 5)
        ['AUGGAGGGGGA', 'AUGGGGAGGGA', 'AGGGGAUGGGA', 'AUGGGGGAGGA', 'AUGGGAGGGGA']
        """
        return [din_s(seq) for _ in range(n)]

    @staticmethod
    def to_fasta(sequences: List[str]) -> str:
        """Combines a list of sequences into a string in FASTA format

        >>> IntaRNApvalue.to_fasta(['AUGGAGGGGGA', 'AUGGGGAGGGA', 'AGGGGAUGGGA', 'AUGGGGGAGGA', 'AUGGGAGGGGA'])
        '>0\\nAUGGAGGGGGA\\n>1\\nAUGGGGAGGGA\\n>2\\nAGGGGAUGGGA\\n>3\\nAUGGGGGAGGA\\n>4\\nAUGGGAGGGGA\\n'
        """
        fasta_str = ''
        n = 0
        for seq in sequences:
            fasta_str += '>{}\n{}\n'.format(n, seq)
            n += 1
        return fasta_str

    """
    def get_score(self, seq_tuple: Tuple[str, str]) -> float:
        Gets the IntaRNA score of a single query/target tuple, query is first element, target second
        res = os.popen('{} -q {} -t {}'.format(self.bin, seq_tuple[0], seq_tuple[1])).read()
        if 'no favorable interaction for target and query' in res:
            return 0.0  # TODO: What do with query/target combinations that don't interact?
        return float(res.split('interaction energy = ')[1].split(' kcal/mol\n')[0])

    def get_scores(self, seq_list: List[Tuple[str, str]]) -> List[float]:
        Gets the IntaRNA score to a given list of query/target tuples, returns a sorted list of scores
        pool = mp.Pool(processes=mp.cpu_count())  # we will do this in parallel and not sequential
        scores = pool.map(self.get_score, seq_list)
        scores.sort()  # Timsort is O(n log n) on average and worst, O(n) on best case
        return scores
    """

    def get_scores(self) -> List[float]:
        """Calculates n IntaRNA scores from random sequences"""
        if self.shuffle_query and self.shuffle_target:
            # TODO: both shuffled
            pass
        else:
            scores = []
            missing = self.n
            while missing:
                if self.shuffle_query:
                    shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing)).encode('utf-8')
                    p = Popen([self.bin, '-q', 'STDIN', '-t', self.target, '--outMode=C', '--outCsvCols=E',
                               '--threads', self.threads], stdout=PIPE, stdin=PIPE)
                else:
                    shuffles = self.to_fasta(self.shuffle_sequence(self.target, missing)).encode('utf-8')
                    p = Popen([self.bin, '-q', self.query, '-t', 'STDIN', '--outMode=C', '--outCsvCols=E',
                               '--threads', self.threads], stdout=PIPE, stdin=PIPE)
                res = p.communicate(input=shuffles)[0].decode().split('\n')  # send shuffles as STDIN
                del res[0], res[-1]  # remove first element aka 'E' and trailing newline element
                scores.extend(res)  # add elements to scores
                missing = self.n - len(scores)
            del scores[self.n-1:-1]  # delete unwanted scores so we have exact amount
            return [float(x) for x in scores]  # return a new lists with all elements as float

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
    import doctest
    doctest.testmod()

    q = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    t = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
        'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    i = IntaRNApvalue()
    i.process_cmd_args(['--query', q, '--target', t, '--amount', '1000', '--shuffle', 'q'])
    t = time.time()
    print(i.get_scores())
    print(time.time() - t)
    # print(i.calculate_pvalue())
    # i.get_dist_function(q, t, 1000, True, True)
    # i.get_pvalue_graph(q, t, 4)
    # print('Integriert: {}'.format(i.calculate_pvalue(q, t, 1000, True, True)))
    # print('Empirisch:  {}'.format(i.calculate_pvalue_empirical(q, t, 1000, True, True)))
