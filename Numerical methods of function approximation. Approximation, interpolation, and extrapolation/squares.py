import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot


x = sp.symbols('x')
X = np.array([0.184, 0.519, 0.854, 1.188, 1.523, 1.858, 2.192, 2.527, 2.862])
Y = np.array([-1.679, -3.056, -4.493, -6.928, -10.524, -15.4, -22.049, -32.38, -44.538])


def get_sympy_subplots(plot: Plot):
    backend = MatplotlibBackend(plot)

    backend.process_series()
    backend.fig.tight_layout()
    return backend.plt


def least_sq(X, Y, k):
    if k >= len(X):
        return
    M = np.zeros((k+1, k+1))
    V = np.zeros(k+1)
    power = 0
    for i in range(k+1):
        for j, powi in zip(range(k+1), range(power, k+1+power)):
            M[i][j] = np.sum(np.power(X, powi))
        power += 1
    M[0][0] = k+1
    for i in range(k+1):
        V[i] = np.sum(np.dot(Y, np.power(X, i)))
    A = np.dot(np.linalg.inv(M), V)
    poltemp = []
    for i in range(k+1):
        poltemp.append(f"({A[i]}*(x**{i}))")
    pol = sp.sympify('+'.join(poltemp))
    Ypol = []
    for i in X:
        Ypol.append(pol.evalf(subs={x: i}))
    sumsqdeltaY = np.sum(np.power(Y-Ypol, 2))
    av_sq_error = np.sqrt(float((sumsqdeltaY)/(len(X)-1)))
    pl = sp.plotting.plot(pol, (x, X[0]-1, X[-1]+1), show=False, c="purple")  # type: ignore
    plt = get_sympy_subplots(pl)
    plt.plot(X, Y, "o", c='red')
    plt.errorbar(X, Ypol, av_sq_error, c='black')
    plt.plot(X, Ypol, ".", c='black')
    plt.show()
    return pol, av_sq_error


plt.scatter(X, Y)
plt.show()
pol, er = least_sq(X, Y, 2)  # type: ignore
print("Aпроксимуюча функція: \n",
      pol)
print("Cередньо-квадратичнa похибкa: \n",
      er)


#print(list(zip(X, Y)))
