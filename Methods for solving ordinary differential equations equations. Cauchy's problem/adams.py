import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
# import pylab as pl

x, y = sp.symbols('x y')


def adams(df, x0, y0, a, b, h):
    Xs = np.arange(a, b+h, h)
    Ys = [y0]
    for i in range(1, 3):
        k1 = df.evalf(subs={x: Xs[i-1], y: Ys[i-1]})
        k2 = df.evalf(subs={x: (Xs[i-1]+(h/2)), y: (Ys[i-1]+(h/2)*k1)})
        k3 = df.evalf(subs={x: (Xs[i-1]+h), y: (Ys[i-1]-h*k1+2*h*k2)})
        Ys.append(Ys[i-1]+(h/6)*(k1+4*k2+k3))

    for i in range(3, len(Xs)):
        k1 = (23*df.evalf(subs={x: Xs[i-1], y: Ys[i-1]}))
        k2 = (16*df.evalf(subs={x: Xs[i-2], y: Ys[i-2]}))
        k3 = (5*df.evalf(subs={x: Xs[i-3], y: Ys[i-3]}))
        Ys.append(Ys[i-1]+(h/12)*(k1-k2+k3))

    print("Метод Адамса третього порядку\n")
    [print(k) for k in list(zip(Xs, Ys))]

    fig, ax = plt.subplots()

    ax.plot(Xs, Ys, 'o')
    ax.plot(Xs, Ys, '-')

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.show()


def main():
    df = 1-sp.sin(x+y)+(0.5*y)/(x+2)
    x0 = 0
    y0 = 0
    a = 0
    b = 1
    h = 0.1
    adams(df, x0, y0, a, b, h)


if __name__ == "__main__":
    main()
