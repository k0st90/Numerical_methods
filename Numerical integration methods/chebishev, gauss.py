import sympy as sp
import numpy as np

x = sp.symbols('x')

Ts = {
    1: [0],
    2: [-0.5773503, 0.5773503],
    3: [-0.7071068, 0, 0.7071068],
    4: [-0.7946545, -0.1875925, 0.1875925, 0.7946545],
    5: [-0.8324975, -0.3745414, 0, 0.3745414, 0.8324975],
    6: [-0.8662468, -0.4225187, -0.2666354, 0.2666354, 0.4225187, 0.8662468],
    7: [-0.8838617, -0.5296568, -0.3239118, 0, 0.3239118, 0.5296568, 0.8838617],
    9: [-0.9115893, -0.6010187, -0.5287618, -0.167906, 0, 0.167906, 0.5287618, 0.6010187, 0.9115893]
}

Ts_for_gauss = {
    1: [0.5],
    2: [-0.5773502, 0.5773502],
    3: [-0.7745967, 0, 0.7745967],
    4: [-0.8611363, -0.3399810, 0.3399810, 0.8611363],
    5: [-0.9061798, -0.5384693, 0, 0.5384693, 0.9061798],
    6: [-0.9324695, -0.6612094, -0.2386192, 0.2386192, 0.6612094, 0.9324695],
    7: [-0.949108, -0.741531, -0.405845, 0, 0.405845, 0.741531, 0.949108],
    8: [-0.960290, -0.796666, -0.525532, -0.183434, 0, 0.183434, 0.525532, 0.796666, 0.960290]
}

Cs = {
    1: [2],
    2: [1, 1],
    3: [0.5555556, 0.8888889, 0.5555556],
    4: [0.3478548, 0.6521452, 0.6521452, 0.3478548],
    5: [0.2369269, 0.4786287, 0.5688889, 0.4786287, 0.2369269],
    6: [0.1713245, 0.3607616, 0.4679139, 0.4679139, 0.3607616, 0.1713245],
    7: [0.129485, 0.279705, 0.38183, 0.417960, 0.38183, 0.279705, 0.129485],
    8: [0.101228, 0.222381, 0.313707, 0.362684, 0.362684, 0.313707, 0.222381, 0.101228]
}


def abs_error_runge_cheb(a, b, n):
    k = n-1
    return (a-b)/((2**k)-1)


def rel_error_runge_cheb(a, b, n):
    return np.abs(abs_error_runge_cheb(a, b, n)/a)*100


def abs_error_runge_gauss(a, b, n):
    k = 2*n-1
    return (a-b)/((2**k)-1)


def rel_error_runge_gauss(a, b, n):
    return np.abs(abs_error_runge_gauss(a, b, n)/a)*100


def cheb(f, left, right, n, T):
    assert (1 <= n <= 9 and n != 8), "n має бути від 1 до 9 і не 8"
    Xs = [((right+left)/2)+((right-left)/2)*t for t in T[n]]
    summ = sum([f.evalf(subs={x: X}) for X in Xs])
    return ((right-left)/n)*summ


def cheb_double(f, left, right, n, T):
    return cheb(f, left, (left+right)/2, n, T) + cheb(f, (left+right)/2, right, n, T)


def gauss(f, left, right, n, T, C):
    assert 1 <= n <= 8, "n має бути від 1 до 8"
    summ = sum([(c*f.evalf(subs={x: ((right+left)/2)+((right-left)/2)*t})) for t, c in zip(T[n], C[n])])
    return ((right-left)/2)*summ


def gauss_double(f, left, right, n, T, C):
    return gauss(f, left, (left+right)/2, n, T, C) + gauss(f, (left+right)/2, right, n, T, C)


def main():
    f = (sp.sqrt(0.3*(x**2)+2.3))/(1.8+sp.sqrt(2*x+1.6))
    f2 = (sp.cos((x**2)+0.6))/(1.2+sp.sin((0.7)*x+0.2))
    left_bound1 = 0.8
    right_bound1 = 1.6
    left_bound2 = 0.5
    right_bound2 = 1.8

    for n in range(1, 10):
        if n != 8:
            print(f"\nЗначення 1 інтегралу за квадратурними формулами Чебишова з h (n={n}):")
            ce = cheb(f, left_bound1, right_bound1, n, Ts)
            print(ce)
            print(f"Значення 1 інтегралу за квадратурними формулами Чебишова з h/2 (n={n}):")
            ced = cheb_double(f, left_bound1, right_bound1, n, Ts)
            print(ced)
            print(f"Абсолютна похибка за Рунге (n={n}):")
            print(abs_error_runge_cheb(ce, ced, n))
            print(f"Відносна похибка (n={n}):")
            print(rel_error_runge_cheb(ce, ced, n), "%")
    print("-"*60)
    for n in range(1, 9):
        print(f"\nЗначення 2 інтегралу за квадратурними формулами Гауса з h (n={n}):")
        ge = gauss(f2, left_bound2, right_bound2, n, Ts_for_gauss, Cs)
        print(ge)
        print(f"Значення 2 інтегралу за квадратурними формулами Гауса з h/2 (n={n}):")
        ged = gauss_double(f2, left_bound2, right_bound2, n, Ts_for_gauss, Cs)
        print(ged)
        print(f"Абсолютна похибка за Рунге (n={n}):")
        print(abs_error_runge_gauss(ged, ge, n))
        print(f"Відносна похибка (n={n}):")
        print(rel_error_runge_gauss(ged, ge, n), "%")




if __name__ == "__main__":
    main()



"""
def simpson(f, left, right, e):
    n = find_n_for_simpson(f, left, right, e)
    h = (right-left)/(n)
    sumf = 0
    prikol = np.linspace(left, right, n+1)
    counter = 0
    for xi in prikol[1:-1]:
        if counter % 2 == 1:
            sumf += 2*f.evalf(subs={x: xi})
        else:
            sumf += 4*f.evalf(subs={x: xi})
        counter += 1

    sumf += f.evalf(subs={x: prikol[0]})
    sumf += f.evalf(subs={x: prikol[-1]})
    return (h/3)*sumf

"""
