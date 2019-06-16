#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

from subprocess import PIPE, Popen, run
import random
import os
import sys
import argparse
from typing import List, Tuple
import time
import numpy as np
from scipy.integrate import quad as integ
from DinuclShuffle import dinucl_shuffle as din_s


class IntaRNApvalue:
    def __init__(self, args=None):
        self.bin = self.find_binary()
        self.query = ''
        self.target = ''
        self.n = 0
        self.shuffle_query = False
        self.shuffle_target = False
        self.threads = ''
        self.process_cmd_args(args)
        self.original_score = self.get_single_score(self.query, self.target)  # TODO: leave this here?
        if not self.original_score:
            print('The query/target combination you specified has no favorable interaction')
            sys.exit(1)

    @staticmethod
    def find_binary() -> str:
        """Tries to find the IntaRNA executable and returns path to binary or exits with error code 1 if not found"""
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

        >>> i = IntaRNApvalue(['-q', 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU',
        ... '-t', 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUU', '--amount', '10', '--shuffle', 'b', '--threads', '3'])
        >>> i.query
        'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
        >>> i.target
        'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUU'
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

    def get_scores(self) -> Tuple[List[float], int]:
        """Calculates n IntaRNA scores from random sequences with given parameters as class variables

        # >>> i = IntaRNApvalue()
        # >>> i.n, i.query, i.target, i.shuffle_query, i.shuffle_target, i.threads = 13, 'ACG', 'GGAU', True, True, '0'
        # >>> i.get_single_score(i.query, i.target)
        """
        scores = []
        missing = self.n
        non_interactions = 0

        while missing > 0:
            if self.shuffle_query and self.shuffle_target:  # shuffle both
                query = self.shuffle_sequence(self.query, 1)[0]  # get a random query
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing)).encode('utf-8')
            elif self.shuffle_query and not self.shuffle_target:  # only shuffle query
                query = 'STDIN'
                target = self.target  # target stays the same
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing)).encode('utf-8')
            else:  # only shuffle target
                query = self.query  # query stays the same
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.target, missing)).encode('utf-8')

            p = Popen([self.bin, '-q', query, '-t', target, '--outMode=C', '--outCsvCols=E', '--threads', self.threads],
                      stdout=PIPE, stdin=PIPE)
            stdout, stderr = p.communicate(input=shuffles)  # send shuffles as STDIN
            stdout = stdout.decode().split('\n')  # convert bytes -> str and split on newline
            del stdout[0], stdout[-1]  # remove first element aka 'E' and trailing newline element
            scores.extend(stdout)  # add elements to scores
            missing = self.n - len(scores)
            non_interactions += missing  # count non-interactions

        # return list with all elements as float and amount of non-interactions
        return [float(x) for x in scores], non_interactions

    def get_single_score(self, query: str, target: str) -> float:
        """Gets an IntaRNA score to a single query/target combination"""
        o = run('{} -q {} -t {} --outMode=C --outCsvCols=E --threads {}'.format(self.bin, query, target, self.threads),
                stdout=PIPE, stdin=PIPE, shell=True).stdout.decode()
        if o.startswith('E'):
            return float(o.split('\n')[1])
        else:
            return 0

    def calculate_pvalue_empirical(self) -> float:
        """Calculates a p-value to a target/query combination empirical with a given amount of shuffle iterations"""
        scores = self.get_scores()[0]
        return [score <= self.original_score for score in scores].count(True) / len(scores)

    def calculate_pvalue(self) -> any:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations"""
        scores = self.get_scores()[0]

        # Try to fit a gaussian distribution
        avg = np.mean(scores)  # average
        var = np.var(scores)  # variance
        # std_dev = np.sqrt(var)  # standard deviation

        def gauss(x):
            return 1.0 / np.sqrt(2 * np.pi * var) * np.exp(-0.5 * ((x - avg) ** 2 / var))
        return integ(gauss, -np.inf, self.original_score)


if __name__ == '__main__':
    i = IntaRNApvalue()

    start = time.time()
    scores = i.get_scores()
    print('We have {} scores and {} non-interactions'.format(len(scores[0]), scores[1]))
    # print(i.calculate_pvalue())
    # print(i.calculate_pvalue_empirical())
    print('This run took: {}'.format(time.time() - start))
    # print(i.calculate_pvalue())
    # i.get_dist_function(q, t, 1000, True, True)
    # i.get_pvalue_graph(q, t, 4)
    # print('Integriert: {}'.format(i.calculate_pvalue(q, t, 1000, True, True)))
    # print('Empirisch:  {}'.format(i.calculate_pvalue_empirical(q, t, 1000, True, True)))
