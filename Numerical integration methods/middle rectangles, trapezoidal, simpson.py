import sympy as sp
import numpy as np
from sympy import solveset, Interval, Max

x = sp.symbols('x')


def find_n_for_midrec(f, lower_bound, upper_bound, e):
    function = sp.diff(f, x, 2).simplify()
    zeros = solveset(function, x, domain=Interval(lower_bound, upper_bound))
    assert zeros.is_FiniteSet
    max_d = Max(sp.Abs(function.subs(x, lower_bound)), sp.Abs(function.subs(x, upper_bound)), *[sp.Abs(function.subs(x, i)) for i in zeros])
    n = sp.Pow((max_d*((upper_bound-lower_bound)**3))/(24*e), 1/2)
    return sp.ceiling(n)


def find_n_for_trapezoidal(f, lower_bound, upper_bound, e):
    function = sp.diff(f, x, 2).simplify()
    zeros = solveset(function, x, domain=Interval(lower_bound, upper_bound))
    assert zeros.is_FiniteSet
    max_d = Max(sp.Abs(function.subs(x, lower_bound)), sp.Abs(function.subs(x, upper_bound)), *[sp.Abs(function.subs(x, i)) for i in zeros])
    n = sp.Pow((max_d*((upper_bound-lower_bound)**3))/(12*e), 1/2)
    return sp.ceiling(n)


def find_n_for_simpson(f, lower_bound, upper_bound, e):
    function = sp.diff(f, x, 4).simplify()
    #print(function)
    zeros = solveset(function, x, domain=Interval(lower_bound, upper_bound))
    assert zeros.is_FiniteSet
    max_d = Max(sp.Abs(function.subs(x, lower_bound)), sp.Abs(function.subs(x, upper_bound)), *[sp.Abs(function.subs(x, i)) for i in zeros])
    n = np.power((max_d*((upper_bound-lower_bound)**5))/(180*e), 1/4)
    if sp.ceiling(n) % 2 == 0:
        return sp.ceiling(n)
    return sp.ceiling(n)+1


def midrec(f, left, right, e):
    n = find_n_for_midrec(f, left, right, e)
    h = (right-left)/n
    sumf = 0
    for xi in np.linspace(left, right, n+1)[:-1]:
        sumf += f.evalf(subs={x: (xi+h/2)})
    return h*sumf


def trapezoidal(f, left, right, e):
    n = find_n_for_trapezoidal(f, left, right, e)
    h = (right-left)/n
    sumf = 0
    prikol = np.linspace(left, right, n+1)
    for xi in prikol[1:-1]:
        sumf += 2*f.evalf(subs={x: xi})
    sumf += f.evalf(subs={x: prikol[0]})
    sumf += f.evalf(subs={x: prikol[-1]})
    return (h/2)*sumf


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


def abs_error(a, b):
    return np.abs(a-b)


def rel_error(a, b):
    return np.abs((a-b)/a)*100


def main():
    f = ((4*x)*(sp.exp(2*x)))/((1+2*x)**2)
    F = sp.exp(2*x)/(1+2*x)
    left_bound = 1
    right_bound = 2
    e1 = 1e-5
    e2 = 1e-6

    true_result = F.evalf(subs={x: right_bound})-F.evalf(subs={x: left_bound})
    print("Точне значення первісної F(x):", true_result)

    midrec_res = midrec(f, left_bound, right_bound, e1)
    trapezoidal_res = trapezoidal(f, left_bound, right_bound, e1)
    simpson_res = simpson(f, left_bound, right_bound, e1)
    print(f"\nЗначення інтегралу за формулою центральних прямокутників (точність {e1}): ")
    print(midrec_res)
    print("Абсолютна похибка: ", abs_error(true_result, midrec_res))
    print("Відносна похибка: ", rel_error(true_result, midrec_res), "%")
    print(f"\nЗначення інтегралу за формулою трапецій (точність {e1}): ")
    print(trapezoidal_res)
    print("Абсолютна похибка: ", abs_error(true_result, trapezoidal_res))
    print("Відносна похибка: ", rel_error(true_result, trapezoidal_res), "%")
    print(f"\nЗначення інтегралу за формулою Сімпсона (точність {e1}): ")
    print(simpson_res)
    print("Абсолютна похибка: ", abs_error(true_result, simpson_res))
    print("Відносна похибка: ", rel_error(true_result, simpson_res), "%")

    midrec_res2 = midrec(f, left_bound, right_bound, e2)
    trapezoidal_res2 = trapezoidal(f, left_bound, right_bound, e2)
    simpson_res2 = simpson(f, left_bound, right_bound, e2)
    print(f"\nЗначення інтегралу за формулою центральних прямокутників (точність {e2}): ")
    print(midrec_res2)
    print("Абсолютна похибка: ", abs_error(true_result, midrec_res2))
    print("Відносна похибка: ", rel_error(true_result, midrec_res2), "%")
    print(f"\nЗначення інтегралу за формулою трапецій (точність {e2}): ")
    print(trapezoidal_res2)
    print("Абсолютна похибка: ", abs_error(true_result, trapezoidal_res2))
    print("Відносна похибка: ", rel_error(true_result, trapezoidal_res2), "%")
    print(f"\nЗначення інтегралу за формулою Сімпсона (точність {e2}): ")
    print(simpson_res2)
    print("Абсолютна похибка: ", abs_error(true_result, simpson_res2))
    print("Відносна похибка: ", rel_error(true_result, simpson_res2), "%")




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
