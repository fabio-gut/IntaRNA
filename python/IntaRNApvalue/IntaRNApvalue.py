#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

from subprocess import PIPE, Popen, run
import random
import os
import argparse
from typing import List, Tuple
import numpy as np
from scipy.integrate import quad as integ
from scipy.stats import norm as gauss
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_l as gum
from IntaRNApvalue.DinuclShuffle import dinucl_shuffle


class IntaRNApvalue:
    def __init__(self, test_args=None):
        self.bin = self.find_binary()
        self.query = ''
        self.target = ''
        self.n = 0
        self.shuffle_query = False
        self.shuffle_target = False
        self.threads = ''
        self.dist = ''
        self.scores_out = False
        self.process_cmd_args(test_args)
        self.original_score = self.get_original_score()
        if not self.original_score and not test_args:  # exit if given seq has no interaction and not test mode
            print('The query/target combination you specified has no favorable interaction')
            exit(1)
        if not test_args:
            self.main()

    def main(self):
        """The main function"""
        scores, non_interactions = self.get_scores()
        if self.scores_out:  # output scores and exits the process
            print('\n'.join(iter([str(x) for x in scores])))
            exit(0)

        pvalue = ''
        if self.dist == 'gauss':
            pvalue = self.calculate_pvalue_gauss(scores)
        elif self.dist == 'none':
            pvalue = self.calculate_pvalue_empirical(scores)
        elif self.dist == 'gumbel':
            pvalue = self.calculate_pvalue_gumbel(scores)
        elif self.dist == 'gev':
            pvalue = self.calculate_pvalue_generalized_ext_val(scores)

        print('pvalue: {}'.format(pvalue))

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
        exit(1)

    def process_cmd_args(self, test_args=None) -> None:
        """Processes all commandline args

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '-n', '10', '-m', 'b', '-d', 'gauss', '-s',
        ... '--threads', '3'])
        >>> i.query
        'AGGAUG'
        >>> i.target
        'UUUAUCGUU'
        >>> i.n
        10
        >>> i.shuffle_query
        True
        >>> i.shuffle_target
        True
        >>> i.threads
        '3'
        >>> i.dist
        'gauss'
        >>> i.scores_out
        True
        """
        parser = argparse.ArgumentParser(description='Calculates p-values to IntaRNA scores.')
        parser.add_argument('-q', '--query', dest='query', type=str, help='Query sequence', required=True)
        parser.add_argument('-t', '--target', dest='target', type=str, help='Target sequence', required=True)
        parser.add_argument('-n', '--scores', dest='n', type=int, required=True,
                            help='How many randomly generated scores are used to calculate the p-value.')
        parser.add_argument('-m', '--shuffle-mode', dest='sm', required=True, choices=['q', 't', 'b'],
                            help='Which sequences are going to be shuffled: both, query only or target only.')
        parser.add_argument('-d', '--distribution', dest='dist', choices=['gev', 'gumbel', 'gauss'], default='gev',
                            help='Which distribution is fitted and used to calculate the pvalue.')
        parser.add_argument('-s', '--scores-out', dest='scores_out', action='store_true',
                            help='All IntaRNA scores used for pvalue calculation are output to STDOUT. '
                                 'There is no pvalue output, so this function is useful for pipeing the scores.')
        # TODO: correct?
        parser.add_argument('--threads', type=str, default='0', help='Sets the amount of threads used for IntaRNA.')
        parser.add_argument('--seed', type=str, default=None,
                            help='Random seed to make sequence generation deterministic.')

        args = parser.parse_args(test_args)

        allowed = ['G', 'A', 'C', 'U']
        args.query = args.query.upper()
        args.target = args.target.upper()
        if False in [n in allowed for n in args.query] or False in [n in allowed for n in args.target]:
            print('A sequence you specified contains illegal characters, allowed: G, A, C, U')
            exit(1)

        self.shuffle_query = True if args.sm in ['b', 'q'] else False
        self.shuffle_target = True if args.sm in ['b', 't'] else False
        for key, value in args.__dict__.items():
            self.__setattr__(key, value)

        random.seed(a=args.seed)

    @staticmethod
    def shuffle_sequence(seq: str, n: int) -> List[str]:
        """Shuffles a sequence n times and returns a list of sequences, duplicate entries are possible

        >>> random.seed('IntaRNA')
        >>> IntaRNApvalue.shuffle_sequence('AGGAUGGGGGA', 5)
        ['AUGGAGGGGGA', 'AUGGGGAGGGA', 'AGGGGAUGGGA', 'AUGGGGGAGGA', 'AUGGGAGGGGA']
        """
        return [dinucl_shuffle(seq) for _ in range(n)]

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
        """Calculates n IntaRNA scores from random sequences with given parameters as class variables"""
        scores = []
        missing = self.n
        non_interactions = 0

        while missing > 0:
            if self.shuffle_query and self.shuffle_target:  # shuffle both
                query = self.shuffle_sequence(self.query, 1)[0]  # get a random query
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing))
            elif self.shuffle_query and not self.shuffle_target:  # only shuffle query
                query = 'STDIN'
                target = self.target  # target stays the same
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing))
            else:  # only shuffle target
                query = self.query  # query stays the same
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.target, missing))

            p = Popen([self.bin, '-q', query, '-t', target, '--outMode=C', '--outCsvCols=E', '--threads', self.threads],
                      stdout=PIPE, stdin=PIPE, universal_newlines=True)
            stdout, stderr = p.communicate(input=shuffles)  # send shuffles as STDIN
            stdout = stdout.split('\n')  # split on newline
            del stdout[0], stdout[-1]  # remove first element aka 'E' and trailing newline element
            scores.extend(stdout)  # add elements to scores
            missing = self.n - len(scores)
            non_interactions += missing  # count non-interactions

        # return list with all elements as float and amount of non-interactions
        return [float(x) for x in scores], non_interactions

    def get_original_score(self) -> float:
        """Gets an IntaRNA score to a single query/target combination"""
        o = run('{} -q {} -t {} --outMode=C --outCsvCols=E --threads {}'.format(self.bin, self.query, self.target,
                                                                                self.threads),
                stdout=PIPE, stdin=PIPE, shell=True).stdout.decode()
        if o.startswith('E') and o != 'E\n':
            return float(o.split('\n')[1])
        else:
            return 0  # no interaction

    def calculate_pvalue_empirical(self, scores: list = None) -> float:
        """Calculates a p-value to a target/query combination empirical with a given amount of shuffle iterations

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-sm', 'b', '--threads', '3'])
        >>> i.original_score = -10.0
        >>> i.calculate_pvalue_empirical([-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.25
        """
        if not scores:
            scores, non_interactions = self.get_scores()
        return [score <= self.original_score for score in scores].count(True) / len(scores)

    def calculate_pvalue_gauss(self, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a gaussian distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-sm', 'b', '--threads', '3'])
        >>> i.original_score = -10.0
        >>> i.calculate_pvalue_gauss([-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.2429106747265256
        """
        loc, scale = gauss.fit(scores)

        def f(x):
            return gauss.pdf(x, loc=loc, scale=scale)
        return integ(f, -np.inf, self.original_score)[0]

    def calculate_pvalue_gumbel(self, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a gumbel distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-sm', 'b', '--threads', '3'])
        >>> i.original_score = -10.0
        >>> i.calculate_pvalue_gumbel([-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.19721934073203196
        """
        loc, scale = gum.fit(scores)

        def f(x):
            return gum.pdf(x, loc=loc, scale=scale)
        return integ(f, -np.inf, self.original_score)[0]

    def calculate_pvalue_generalized_ext_val(self, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a generalized extreme value distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-sm', 'b', '--threads', '3'])
        >>> i.original_score = -10.0
        >>> i.calculate_pvalue_generalized_ext_val([-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.17611816922560236
        """
        shape, loc, scale = gev.fit(scores)

        def f(x):
            return gev.pdf(x, shape, loc=loc, scale=scale)
        return integ(f, -np.inf, self.original_score)[0]


if __name__ == '__main__':
    i = IntaRNApvalue()
