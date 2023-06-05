import sympy as sp
import numpy as np


X = [0.015, 0.681, 1.342, 2.118, 2.671]
Y = [-2.417, -3.819, -0.642, 0.848, 2.815]

x = sp.symbols('x')


def lag(X, Y):
    summ = []
    for j in range(len(X)):
        temp = []
        for m in range(len(X)):
            if m != j:
                temp.append(f"((x-{X[m]})/({X[j]}-{X[m]}))")
        summ.append(f"{Y[j]}*{'*'.join(temp)}")
    res = (sp.expand(sp.sympify("+".join(summ))))
    sp.plotting.plot(res, xlim=(-2, 6), ylim=(-8, 12), autoscale=True)  # type: ignore
    return res


def linear_spline(X, Y):
    outlist = []
    for i in range(len(X)-1):
        outlist.append(sp.sympify(f"{Y[i]}+"+f"(({Y[i+1]}-{Y[i]})/({X[i+1]}-{X[i]}))*(x-{X[i]})"))
    p1 = sp.plotting.plot(outlist[0],
                          (x, X[0], X[1]),
                          show=False,
                          autoscale=False,
                          backend='matplotlib',
                          markers=[{'args': [X, Y, 'ro']}])  # type: ignore
    count = 1
    for eq in outlist[1:]:
        p1.append(sp.plotting.plot(eq, (x, X[count], X[count+1]), show=False)[0])  # type: ignore
        count += 1
    p1.show()
    return outlist


def qudratic_spline(X, Y):
    outlist = []
    m = (len(X)-1)*3
    M = np.zeros((m, m))
    b = np.zeros(m)
    counter = 0
    for i, j in zip(range(0, (len(X)-1)*2, 2), range(0, m, 3)):
        M[i][j] = np.square(X[counter])
        M[i][j+1] = X[counter]
        M[i][j+2] = 1
        M[i+1][j] = np.square(X[counter+1])
        M[i+1][j+1] = X[counter+1]
        M[i+1][j+2] = 1
        counter += 1

    counterd = 1
    for i, j in zip(range((len(X)-1)*2, m-1, 1), range(0, m, 3)):
        M[i][j] = 2*X[counterd]
        M[i][j+1] = 1
        M[i][j+3] = -2*X[counterd]
        M[i][j+4] = -1
        counterd += 1
    M[(len(X)-1)*3-1][0] = 2

    b[0] = Y[0]
    b[len(Y)*2-3] = Y[-1]

    for k, y in zip(range(1, len(Y)*2-3, 2), Y[1:]):
        b[k] = y
        b[k+1] = y

    abc = np.dot(np.linalg.inv(M), b)
    for a in range(0, len(abc), 3):
        outlist.append(sp.sympify(f"({abc[a]}*x**2)+({abc[a+1]}*x)+({abc[a+2]})"))

    p1 = sp.plotting.plot(outlist[0], (x, X[0], X[1]),
                          axis_center=(0, 0),
                          margin=1,
                          show=False,
                          autoscale=True,
                          backend='matplotlib',
                          markers=[{'args': [X, Y, 'ro']}]) # type: ignore
    count = 1
    for eq in outlist[1:]:
        p1.append(sp.plotting.plot(eq, (x, X[count], X[count+1]), show=False)[0])  # type: ignore
        count += 1
    p1.show()
    return outlist

def eval_sp(X, list):
    target_number = X[0] + X[1]
    num_of_sp = 0
    for i in range(len(X) - 1):
        if X[i] < target_number < X[i + 1]:
            num_of_sp = i

    res = list[num_of_sp].evalf(subs={x: target_number})
    return res


def main():

    lagr = lag(X, Y)
    print("Інтерполяційний поліном Лагранжа: ")
    print(lagr)
    lines = linear_spline(X, Y)
    print("Лінійний сплайн: ")
    [print(s) for s in lines]
    qudros = qudratic_spline(X, Y)
    print("Квадтратичний сплайн: ")
    [print(s) for s in qudros]
    print("S1(x1+x2): ", eval_sp(X, lines))
    print("S2(x1+x2): ", eval_sp(X, qudros))
    print("L(x1+x2): ", lagr.evalf(subs={x: X[0] + X[1]}))


if __name__ == "__main__":
    main()
