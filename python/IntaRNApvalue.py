#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

import os
import sys
import argparse
import subprocess
from subprocess import PIPE
from DinuclShuffle import dinucl_shuffle as din_s
import time
from matplotlib import pyplot as plt


class IntaRNApvalue:
    def __init__(self):
        self.bin = self.find_binary()

    @staticmethod
    def find_binary() -> str:
        """Tries to find the IntaRNA executable and returns path to binary or exits with error code 1 if not found"""
        to_search = os.environ['PATH'].split(':')  # get path env variable
        to_search.append(os.path.abspath(os.path.join(os.curdir, '../..')))  # add 2 dirs above this script

        for path in to_search:  # search in each directory including subdirectories
            for dir_path, dir_name, file_names in os.walk(path):
                if 'IntaRNA' in file_names:
                    return os.path.join(dir_path, 'IntaRNA')
        else:
            print('Error: Cannot find IntaRNA binary executable, please add it to your PATH')
            sys.exit(1)

    def process_cmd_args(self):
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
        self.calculate_pvalue(args.query, args.target, args.amount, args.sq, args.st)

    def get_score(self, query: str, target: str) -> float:
        """Gets the IntaRNA score of a single target/query combination"""
        res = os.popen('{} -q {} -t {}'.format(self.bin, query, target)).read()
        # res = subprocess.run([self.bin, '-q{}'.format(query), '-t{}'.format(target)], stdout=PIPE, stderr=PIPE)
        # TODO: For some reason both subprocess.run() and os.popen() are both slow on my system
        # stdout = res.stdout.decode('utf-8')
        # if res.returncode:  # if process ended with errcode != 0
        #     print(res.stderr.decode('utf-8'))
        #     sys.exit(res.returncode)  # exit with same error code
        if 'no favorable interaction for target and query' in res:
            return 0.0  # TODO: What do with query/target combinations that don't interact?
        return float(res.split('interaction energy = ')[1].split(' kcal/mol\n')[0])

    def get_scores(self, seq_list: list) -> list:
        """Gets the IntaRNA score to a given list of query/target tuples, returns a sorted list of scores"""
        scores = []
        for query, target in seq_list:
            scores.append(self.get_score(query=query, target=target))
        scores.sort()  # Timsort is O(n log n) on average and worst, O(n) on best case
        return scores

    @staticmethod
    def shuffle_sequence(query: str, target: str, n: int, shuffle_query: bool, shuffle_target: bool) -> list:
        """Shuffles a query/target pair n times and returns an array of sequence tuples, by default both are shuffled.
            The returned array always has length n+1, duplicate entries are possible"""
        shuffles = [(query, target)]  # TODO: add original seq???
        for _ in range(n):
            shuffles.append((din_s(query) if shuffle_query else query, din_s(target) if shuffle_target else target))
        return shuffles

    def calculate_pvalue(self, query: str, target: str, n: int, shuffle_query: bool, shuffle_target: bool) -> float:
        """Calculates a p-value to a target/query combination with a given amount of shuffle iterations"""
        start = time.time()
        original_score = self.get_score(query, target)
        shuffles = self.shuffle_sequence(query, target, n, shuffle_query, shuffle_target)
        shuffle_start = time.time()
        scores = self.get_scores(shuffles)

        total = time.strftime("%H:%M:%S.%f", time.gmtime(time.time() - start))
        shuffle = time.strftime("%H:%M:%S.%f", time.gmtime(shuffle_start - start))
        calculate = time.strftime("%H:%M:%S.%f", time.gmtime(time.time() - shuffle_start))
        print('This run took {} ({} to shuffle and {} to calculate scores)'.format(
            total, shuffle, calculate
        ))
        # Do it empirical for now, TODO: fit curve
        return [score <= original_score for score in scores].count(True) / len(scores)

    def get_pvalue_graph(self, query: str, target: str, max_exp: int = 4):
        for n in range(1, max_exp):
            plt.plot(10**n, self.calculate_pvalue(query, target, 10**n, True, False), 'rx')
            plt.plot(10**n, self.calculate_pvalue(query, target, 10**n, False, True), 'bx')
            plt.plot(10**n, self.calculate_pvalue(query, target, 10**n, True, True), 'gx')
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
    i.get_pvalue_graph(q, t, 4)
