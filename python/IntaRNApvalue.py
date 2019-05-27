#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

import os
import sys
import argparse
import subprocess
from subprocess import PIPE


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
            print('Error: Cannot find IntaRNA binary executable')
            sys.exit(1)

    def process_cmd_args(self):
        """Processes all commandline args"""
        parser = argparse.ArgumentParser(description='Calculates p-values to IntaRNA scores')
        parser.add_argument(['--query', '-q'], type=str, required=True, help='Query sequence')
        parser.add_argument(['--target', '-t'], type=str, required=True, help='Target sequence')
        parser.add_argument(['--query-shuffles', '-qn'], type=int, default=10000, dest='qn',
                            help='Amount of shuffles to the query sequence')
        parser.add_argument(['--target-shuffles', '-tn'], type=int, default=10000, dest='tn',
                            help='Amount of shuffles to the target sequence')
        args = parser.parse_args()
        self.calculate_pvalue(query=args.query, target=args.target, query_n=args.qn, target_n=args.tn)

    def get_score(self, query: str, target: str) -> float:
        """Gets the IntaRNA score of a given target/query combination"""
        res = subprocess.run([self.bin, '-q{}'.format(query), '-t{}'.format(target)], stdout=PIPE, stderr=PIPE)
        if res.returncode:  # if process ended with errcode != 0
            print(res.stdout.decode('utf-8'))
            sys.exit(res.returncode)
        return float(res.stdout.decode('utf-8').split('interaction energy = ')[1].split(' kcal/mol\n')[0])

    def shuffle_sequence(self, sequence: str, n: int) -> list:
        """Shuffles a sequence n times and returns an array of sequences
        TODO: remove duplicates?"""
        pass

    def calculate_pvalue(self, query: str, target: str, query_n: int, target_n: int) -> float:
        """Calculates a p-value to a target/query combination with a given amount of shuffle iterations"""
        pass


if __name__ == '__main__':
    query = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    target = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
             'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    i = IntaRNApvalue()
    i.process_cmd_args()
    # print(i.get_score(query=query, target=target))
