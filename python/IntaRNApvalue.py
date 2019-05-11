# Copyright 2019
# Author: Fabio Gutmann <fabio.gutmann@jupiter.uni-freiburg.de>

import os


def find_binary():
    """Tries to find the IntaRNA executable, returns None if not found"""

    to_check = [os.path.abspath(os.path.join(os.curdir, '../..')), '/usr/local/bin/', '/usr/bin/']
    for path in to_check:
        for dir_path, dir_name, file_names in os.walk(path):
            if 'IntaRNA' in file_names:
                return os.path.join(dir_path, 'IntaRNA')
    else:
        return None


if __name__ == '__main__':
    query = 'AGGAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGCUUUUUU'
    target = 'UUUAAAUUAAAAAAUCAUAGAAAAAGUAUCGUUUGAUACUUGUGAUUAUACUCAGUUAUACAGUAUCUUAAGGUGUUAUUAAUAGUGGUG' \
             'AGGAGAAUUUAUGAAGCUUUUCAAAAGCUUGCUUGUGGCACCUGCAACUCUUGGUCUUUUAGCACCAAUGACCGCUACUGCUAAU'

    intarna_bin = find_binary()
    if not intarna_bin:
        print('Error: Cannot find IntaRNA binary executable')
        os._exit(1)

    test = os.system('IntaRNA -q {0} -t {1}'.format(query, target))

