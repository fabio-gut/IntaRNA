# IntaRNApvalue
A tool for calculating pvalues to IntaRNA scores.

### How it works:
This tool shuffles one or more of the given sequences (depending on shuffle mode) while preserving mono- and dinucleotide frequences.
Thus the resulting sequences are very similar with similar properties to the original one.
It then calculates the IntaRNA interaction scores for these newly generated sequences and uses them to approximate the pvalue of the original target/query combination.
This pvalue is a score for the likelihood that a sequence with a better interaction score than the original ones can randomly occur.
It can thus be used as a measurement for how good an interaction between two sequences is.

## Dependencies:
- IntaRNA
##### And if you want to run it with python or compile it yourself:
- Python 3
- numpy
- scipy

## Installation:
You can either run this tool directly with python, install it with setuptools or run it as a compiled binary.
##### Run with python directly:
You don't have to install it, if you just plan on running it from a certain directory.
Simply copy the "intarnapvalue" folder where you want to run it from.
You have to get the dependencies yourself in this case.

To use the tool from anywhere however, you need to install it as a python module.
All dependencies (except IntaRNA) will be installed automatically.
This way it can also be installed in a virtual environment.
```sh
python setup.py install
```

##### Run as binary:
Go into the bin directory. You can find pre-compiled builds for linux and windows for x64 directly.
You don't need to install any dependencies for running these.
If you want to compile it yourself, you need to get the dependencies and PyInstaller.
Simply run the build.py with the python binary that has the dependencies and PyInstaller installed.
You will find your binary in the build folder.

## Usage:
Run it as a module without installing (you need to be in the parent dir of intarnapvalue):
```sh
python3 -m intarnapvalue <arguments>
```

Run it as a module from anywhere with python (need to run setup.py first, see [Installation](#installation)):
```shell
python3 -m intarnapvalue <arguments>
```
Or use as a compiled binary:
```console
./IntaRNApvalue <arguments>
```

Or import it and use it from within other python code:
```python
from intarnapvalue.intarna_pvalue import IntaRNApvalue
IntaRNApvalue(['--flag1', 'arg1'])
```

## Arguments:

| Flag               | Value                | Default | Description          |
| ------------------ |:-------------------- | :------ | -------------------- |
| -h, --help         | None                 | None    | Gives detailed command help.  |
| -q, --query        | Sequence             | None    | The query sequence.   |
| -t, --target       | Sequence             | None    | The target sequence.  |
| -s, --scores       | Integer              | None    | How many scores are used from randomly permuted sequences for pvalue calculation. |
| -m, --shuffle-mode | {q, t, b}            | None    | Which sequence will be shuffled: Query, Target or both. |
| -d, --distribution | {gev, gumbel, gauss} | gev     | The distribution used for pvalue calculation: Generalized Extreme Value Distribution, Gumbel Distribution or Gauss. |
| -o, --output       | {pvalue, scores}     | pvalue  | If set to pvalue, will only output pvalue. If set to scores, will output every score from randomly generated sequences, but no pvalue. |
| --threads          | 0 - {max threads}    | 0       | How many threads IntaRNA uses for score calculation. If set to 0 it will use all available threads. |
| --seed             | any                  | None    | The seed used for generating random sequences. |