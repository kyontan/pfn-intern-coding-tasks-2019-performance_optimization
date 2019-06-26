# PFN intern coding tasks (2019, Performance Optimization)

## Environment

- Python 3.7.2
- ply 3.11
- numpy 1.16.2
- scipy 1.2.1

## How to run

- `python opt.py [N] | python interpreter.py [N]`
  - run FFT
- `python opt.py [N] | python cost.py`
  - show calcuration cost
- `python test_opt.py`
  - run unit test

## Result

I implemented Cooley-Tukey FFT and Rader's algorithm.

- Sum of cost from N = 2 to 256 (calculated by `cost.py`)
  - `naive.py` (original): 16875645
  - `opt.py`: 3564224 (21.1%)
- max error
  - 8.83390668259562e−13 (when N=212, initializer=sin)

## License

Unknown, beucase this repository includes source code base on https://github.com/pfnet/intern-coding-tasks

I don't encourage you to include these source codes into other program.
