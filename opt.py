#!/usr/bin/env python3

import argparse

def f(m, n):
    if m == 0:
      return '(1)'
    return 'f({}, {})'.format(m, n)

def e(m):
    return 'e_{}'.format(m)

# return factors of n (including duplicates)
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

# emit optimized version of `{yexp} = {xexp} * f(m, n)`
def emit_y_add_x_prod_f(yexp, xexp, m, n, temp='t'):
    if m % n == 0: # f(m, n) = 1
        print('{} = {} + {}'.format(yexp, yexp, xexp))
    elif (m % n) % n == n/2: # f(m, n) = -1
        print('{} = {} - {}'.format(yexp, yexp, xexp))
    elif (m % n) % n == n/4: # f(m, n) = i
        print('{} = {} * (j)'.format(temp, xexp))
        print('{} = {} + {}'.format(yexp, yexp, temp))
    elif (m % n) % n == n*3/4: # f(m, n) = -i
        print('{} = {} * (-j)'.format(temp, xexp))
        print('{} = {} + {}'.format(yexp, yexp, temp))
    else:
        print('{} = {} * f({}, {})'.format(temp, xexp, m, n))
        print('{} = {} + {}'.format(yexp, yexp, temp))

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

    # Cooley-Tukey's algorithm (1): reorder input array by factors
    prod_factor = 1
    for f in factors[:-1]:
        swap_xi_temp = [None] * N
        m = N // prod_factor
        q = m // f
        for i in range(prod_factor):
            for j in range(f):
                for k in range(q):
                    swap_xi_temp[i*m + q*j + k] = swap_xi_from[i*m + f*k + j]

        swap_xi_from = swap_xi_temp

        prod_factor *= f

    for i in range(N):
        print('y[{}] = x[{}]'.format(i, swap_xi_from[i]))

    step = 0

    n = N // factors[-1]
    m = N // n
    for i in range(n):
        # DFT n sub-arrays: [0..m-1], [m..2m-1], ..., [m(n-1)..mn-1]
        for j in range(m):
            yexp = 'y_{}_{}'.format(step, i*m+j)
            for k in range(m):
                xexp = 'y[{}]'.format(i*m+k) # use y[] as input

                t = j * k

                if k == 0:
                    print('{} = {}'.format(yexp, xexp))
                else:
                    emit_y_add_x_prod_f(yexp, xexp, t, m)

    # Cooley-Tukey's algorithm (2)
    n = N // factors[-1]
    for f in reversed(factors[:-1]):
        step += 1
        m = N // n
        n = n // f
        q = N // n

        for i in range(n):
            for j in range(f):
                for k in range(m):
                    yexp = 'y_{}_{}'.format(step, i*f*m + j*m+k)
                    for l in range(f):
                        lastyexp = 'y_{}_{}'.format(step-1, i*f*m + l*m+k)

                        t = l*(j*m+k)

                        if l == 0:
                            print('{} = {}'.format(yexp, lastyexp))
                        else:
                            emit_y_add_x_prod_f(yexp, lastyexp, t, q)

    for i in range(N):
        print("y[{}] = y_{}_{}".format(i, step, i))

if __name__ == '__main__':
    main()
