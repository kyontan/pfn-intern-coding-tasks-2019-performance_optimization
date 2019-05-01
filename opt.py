#!/usr/bin/env python3

import argparse

def f(m, n):
    if m == 0:
      return '(1)'
    return 'f({}, {})'.format(m, n)

def e(m):
    return 'e_{}'.format(m)

def prime_factorize(n):
    factors = []
    prod_factor = 1
    while prod_factor != n:
        found = False
        for i in range(2, int(n ** 0.5) + 2):
            if (n // prod_factor) % i == 0:
                factors.append(i)
                prod_factor *= i
                found = True
                break
        if not found:
            factors.append(n // prod_factor)
            break

    p = 1
    for f in factors:
        p *= f
    assert(p == n)

    return factors

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='dim', type=int)
    args = parser.parse_args()

    N = args.dim

    if N == 1:
        print("y[0] = x[0]")
        return

    factors = prime_factorize(N)

    swap_xi_from = [x for x in range(N)]

    prod_factor = 1
    for f in factors[:-1]:
        swap_xi_temp = [None] * N
        m = N // prod_factor
        q = m // f
        for i in range(prod_factor):
            for j in range(f):
                for k in range(q):
                    swap_xi_temp[i*m + q*j + k] = swap_xi_from[i*m + f*k + j]

        for i in range(N):
            swap_xi_from = swap_xi_temp

        prod_factor *= f

    # for i in range(N):
    #     print("y[{}] = (0)".format(i))

    step = 0

    # DFT N // factors[-1] sub-arrays
    n = N // factors[-1]
    m = N // n
    for i in range(n):
        for j in range(m):
            yexp = 'y_{}_{}'.format(step, i*m+j)
            for k in range(m):
                xexp = 'x[{}]'.format(swap_xi_from[i*m+k])

                if k == 0:
                    print('{} = {}'.format(yexp, xexp))
                elif (((j*k)%m)*n / N) == 0.5:
                    print('{} = {} - {}'.format(yexp, yexp, xexp))
                else:
                    print('t = {} * f({}, {})'.format(xexp, ((j*k)%m)*(n), N))
                    print('{} = {} + t'.format(yexp, yexp))

    # for i in range(N):
    #     print("x_{} = y[{}]".format(i, i))

    n = prod_factor
    for f in reversed(factors[:-1]):
        step += 1
        m = N // n
        q = n // f
        r = N // q
        for i in range(N):
            print("y_{}_{} = (0)".format(step, i))

        for i in range(q):
            for j in range(f):
                for k in range(m):
                    yexp = 'y_{}_{}'.format(step, i*f*m+(j*m+k))
                    for l in range(f):
                        lastyexp = 'y_{}_{}'.format(step-1, i*f*m + l*m+k)

                        if l == 0:
                            print('{} = {}'.format(yexp, lastyexp))
                        elif ((l*(j*m+k)%r)*q / N) == 0.5:
                            print('{} = {} - {}'.format(yexp, yexp, lastyexp))
                        else:
                            print('t = {} * f({}, {})'.format(lastyexp, (l*(j*m+k)%r)*q, N))
                            print('{} = {} + t'.format(yexp, yexp))

        # print(m, q, r)
        n = q
    for i in range(N):
        print("y[{}] = y_{}_{}".format(i, step, i))

if __name__ == '__main__':
    main()
