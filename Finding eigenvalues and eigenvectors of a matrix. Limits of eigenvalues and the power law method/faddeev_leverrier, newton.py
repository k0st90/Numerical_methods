import numpy as np
import sympy as sp
x = sp.Symbol("x")


def faddeev_leverrier(M):
    M = np.array(M)
    n = M.shape[0]
    assert M.shape[1] == n, 'Array must be square!'

    a = np.array([1.])
    Ak = np.array(M)
    for k in range(1, n + 1):
        ak = -Ak.trace() / k
        a = np.append(a, ak)
        Ak += np.diag(np.repeat(ak, n))
        Ak = np.dot(M, Ak)
    return a


def newton(f, df, x0, tol):
    global x
    if abs(f.evalf(30, subs={x: x0})/df.evalf(30, subs={x: x0})) < tol:
        return x0
    else:
        return newton(f, df, x0 - f.evalf(30, subs={x: x0})/df.evalf(30, subs={x: x0}), tol)


def formp(v):

    global x

    f = sp.sympify("+".join([f"(x**{np.abs(i-np.abs(v.shape[0]-1))})*{a}" for i, a in enumerate(v)]))
    print(f)
    df = sp.diff(f, x)
    return f, df


def solvefln(M):
    v = faddeev_leverrier(M)
    f, df = formp(v)
    sp.plotting.plot(f)  # type: ignore
    listin = []
    for q in range(v.shape[0]-1):
        listin.append(float(input(f"{q+1} Enter: ")))

    outlist = []
    for a in listin:
        out = newton(f, df, a, 0.001)
        outlist.append(out)
    return outlist


def main():
    A = np.array([[1.7, 0.8, 0.9],  # type: ignore
                  [0.8, 0.7, 0.3],
                  [0.9, 0.3, 1.7]])

    print(solvefln(A))


if __name__ == '__main__':
    main()
