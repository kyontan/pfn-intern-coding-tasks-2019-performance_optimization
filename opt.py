#!/usr/bin/env python3

import argparse
from math import log2, ceil

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

def prime_factorize_unique(n):
    return sorted(list(set(prime_factorize(n))))

def prime_primitive_root(p):
    if p == 2:
        return 2
    for i in range(2, p-1):
        for j in prime_factorize_unique(p-1):
            if pow(i, (p-1)//j, p) == 1:
                break
        else:
            return i
    raise "failed to find prime primitive root"

# emit optimized version of `{yexp} = {xexp} * f(m, n)`
def emit_y_add_x_prod_f(yexp, xexp, m, n, temp='t'):
    if m % n == 0: # f(m, n) = 1
        print('{} = {} + {}'.format(yexp, yexp, xexp))
    elif (m % n) % n == n/2: # f(m, n) = -1
        print('{} = {} - {}'.format(yexp, yexp, xexp))
    elif (m % n) % n == n/4: # f(m, n) = i
        print('{} = {} * (-j)'.format(temp, xexp))
        print('{} = {} + {}'.format(yexp, yexp, temp))
    elif (m % n) % n == n*3/4: # f(m, n) = -i
        print('{} = {} * (j)'.format(temp, xexp))
        print('{} = {} + {}'.format(yexp, yexp, temp))
    else:
        print('{} = {} * f({}, {})'.format(temp, xexp, m, n))
        print('{} = {} + {}'.format(yexp, yexp, temp))

def dft(src, dst, size, n=1, sign=1):
    m = size // n
    for i in range(n):
        # DFT n sub-arrays: [0..m-1], [m..2m-1], ..., [m(n-1)..mn-1]
        for j in range(m):
            yexp = '{}_{}'.format(dst, i*m+j)
            for k in range(m):
                lastyexp = '{}_{}'.format(src, i*m+k)

                t = j * k

                if k == 0:
                    print('{} = {}'.format(yexp, lastyexp))
                else:
                    emit_y_add_x_prod_f(yexp, lastyexp, sign * t, m)

def idft(src, dst, size, n=1):
    dft(src, dst, size, n, sign=-1)
    m = size // n
    for i in range(n):
        for j in range(m):
            print("{}_{} = {}_{} * ({})".format(dst, i*m + j, dst, i*m + j, 1.0/m))

def fft(src, dst, size, sign=1):
    # print("fft({}, {}, size={}, sign={})".format(src, dst, size, sign))
    n = 1
    step = 0
    temp = "fft_temp"
    for i in range(size):
        print("{}_{}_{} = {}_{}".format(temp, step, i, src, i))
    step += 1

    while n < size:
        m = size // n
        q = m >> 1
        for i in range(n):
            for k in range(q):
                print("{}_{}_{} = {}_{}_{}".format(temp, step, i*m + k,     temp, step - 1, i*m + 2*k))
                print("{}_{}_{} = {}_{}_{}".format(temp, step, i*m + q + k, temp, step - 1, i*m + 2*k + 1))
        n <<= 1
        step += 1

    while 1 < n:
        m = size // n
        q = n >> 1
        r = size // q

        for i in range(q):
            for k in range(m):
                texp = "{}_{}_{}".format(temp, step, 2*m*i + k)
                texp2 = "{}_{}_{}".format(temp, step, 2*m*i + m + k)
                lasttexp = "{}_{}_{}".format(temp, step - 1, 2*m*i + k)
                lasttexp2 = "{}_{}_{}".format(temp, step - 1, 2*m*i + m + k)
                print("{} = {} * f({}, {})".format(texp, lasttexp2, sign * (k % r) * q, size))
                print("{} = {} + {}".format(texp, texp, lasttexp))
                print("{} = {} * f({}, {})".format(texp2, lasttexp2, sign * ((m + k) % r) * q, size))
                print("{} = {} + {}".format(texp2, texp2, lasttexp))
        step += 1
        n >>= 1

    for i in range(size):
        print("{}_{} = {}_{}_{}".format(dst, i, temp, step - 1, i))

def ifft(src, dst, size):
    fft(src, dst, size, sign=-1)
    for i in range(size):
        print("{}_{} = {}_{} * ({})".format(dst, i, dst, i, 1.0/size))


# Rader's algorithm
def dft_prime_opt(src, dst, size, n=1):
    p = size // n
    g = prime_primitive_root(p)
    ig = pow(g, p - 2, p) # inverse mod of g

    # p-1 <= l && l is power of 2
    l = p - 1
    if not log2(l) % 1 == 0:
        l = pow(2, ceil(log2(p-1)) + 1)
    padlen = l - (p - 1)

    temp = "dft_prime_opt_temp"
    fsrc = '{}_f'.format(temp)
    fdst = '{}_fft_f'.format(temp)
    for i in range(l):
        print("{}_{} = f({}, {})".format(fsrc, i, pow(ig, i%(p-1), p), p))
        # print("{}_{} = f({}, {})".format(fsrc, i, pow(ig, i, p), p))
    fft(fsrc, fdst, size=l)

    for i in range(n):
        step = 1
        for j in range(l):
            texp = "{}_{}_{}".format(temp, step, j)
            if j == 0:
                print("{} = {}_{}".format(texp, src, i*p + 1))
            elif j <= padlen:
                print("{} = (0)".format(texp))
            else:
                print("{} = {}_{}".format(texp, src, i*p + pow(g, (j-padlen) % (p-1),p)))
                # print("{}_{} = {}_{}".format(temp1, j, src, i*p + pow(g, j-padlen, p)))
        step += 1
        fft_src = '{}_{}'.format(temp, step - 1)
        fft_dst = '{}_{}'.format(temp, step)
        fft(fft_src, fft_dst, size=l)

        for j in range(l):
            texp = "{}_{}_{}".format(temp, step, j)
            print("{} = {} * {}_{}".format(texp, texp, fdst, j))

        step += 1
        ifft_src = fft_dst
        ifft_dst = '{}_{}'.format(temp, step)
        ifft(ifft_src, ifft_dst, size=l)

        step += 1
        for j in range(l):
            print("{}_{}_{} = {}_{}_{} + {}_{}".format(temp, step, j, temp, step-1, j, src, i*p))

        # DC (dst[i*p])
        print("{}_{} = {}_{}".format(dst, i*p, src, i*p))
        for j in range(1,p):
            print("{}_{} = {}_{} + {}_{}".format(dst, i*p, dst, i*p, src, i*p + j))

        # not DC (dst[i*p + j])
        for j in range(p-1):
            print("{}_{} = {}_{}_{}".format(dst, i*p + pow(ig, j, p), temp, step, j))

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

    step = 0

    for i in range(N):
        print('y_{}_{} = x[{}]'.format(step, i, swap_xi_from[i]))
    step += 1

    n = N // factors[-1]

    src = 'y_{}'.format(step-1)
    dst = 'y_{}'.format(step)
    if 4 < factors[-1]:
        dft_prime_opt(src, dst, N, n)
    else:
        dft(src, dst, N, n)

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
