import sympy as sp
import numpy as np

x = sp.symbols('x')

Hs = {
    1: [1/2, 1/2],
    2: [1/6, 4/6, 1/6],
    3: [1/8, 3/8, 3/8, 1/8],
    4: [7/90, 32/90, 12/90, 32/90, 7/90],
    5: [19/288, 75/288, 50/288, 50/288, 75/288, 19/288],
    6: [41/840, 216/840, 27/840, 272/840, 27/840, 216/840, 41/840],
    7: [751/17280, 3577/17280, 1323/17280, 2889/17280, 2889/17280, 1323/17280, 3577/17280, 751/17280],
    8: [989/28350, 5888/28350, -928/28350, 10496/28350, -4540/28350, 10496/28350, -928/28350, 5888/28350, 989/28350]
}


def abs_error_runge(a, b, n):
    lk = [1, 3, 3, 5, 5, 7, 7, 9]
    k = lk[n-1]
    return (a-b)/((2**k)-1)


def rel_error_runge(a, b, n):
    return np.abs(abs_error_runge(a, b, n)/a) * 100


def newt_c(f, left, right, n, H):
    assert 1 <= n <= 8, "n має бути від 1 до 8"
    Xs = np.linspace(left, right, n+1)
    summ = sum([h*f.evalf(subs={x: X}) for h, X in zip(H[n], Xs)])
    return (right-left)*summ


def newt_c_double(f, left, right, n, H):
    return newt_c(f, left, (left+right)/2, n, H) + newt_c(f, (left+right)/2, right, n, H)


def main():
    f = (1+0.9*(x**2))/(1.3+sp.sqrt(0.5*(x**2)+1))
    left_bound = 0.9
    right_bound = 2.34
    for n in range(1, 9):
        print(f"\nЗначення інтегралу за квадратурними формулами Ньютона-Котеса з h (n={n}):")
        ne = newt_c(f, left_bound, right_bound, n, Hs)
        print(ne)
        print(f"Значення інтегралу за квадратурними формулами Ньютона-Котеса з h/2 (n={n}):")
        ned = newt_c_double(f, left_bound, right_bound, n, Hs)
        print(ned)
        print(f"Абсолютна похибка за Рунге (n={n}):")
        print(abs_error_runge(ned, ne, n))
        print(f"Відносна похибка (n={n}):")
        print(rel_error_runge(ned, ne, n), "%")





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
