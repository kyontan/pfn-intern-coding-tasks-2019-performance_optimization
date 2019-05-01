#!/usr/bin/env python3

import argparse

def f(m, n):
    if m == 0:
      return '(1)'
    return 'f({}, {})'.format(m, n)

def e(m):
    return 'e_{}'.format(m)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='dim', type=int)
    args = parser.parse_args()

    N = args.dim

    already_defined = [False] * (N * N)

    for l in range(N):
        for k in range(N):
            if already_defined[k*l]:
                pass
            else:
                print('e_{} = {}'.format(k * l, f(k * l, N)))
                already_defined[k*l] = True

    for l in range(N):
        print('t = (0)')
        for k in reversed(range(N)):
            print('t = t + x[{}]'.format(k))
            if k != 0 and l != 0:
                print('t = t * {}'.format(e(l)))

        print('y[{}] = t'.format(l, l))

if __name__ == '__main__':
    main()
