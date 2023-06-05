import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
# import pylab as pl

x, y = sp.symbols('x y')


def runge_kutta_4(df, x0, y0, a, b, h):
    Xs = np.arange(a, b+h, h)
    Ys = [y0]
    # print(Xs[0], Ys[0])
    for i in range(1, len(Xs)):
        k1 = df.evalf(subs={x: Xs[i-1], y: Ys[i-1]})
        k2 = df.evalf(subs={x: (Xs[i-1]+(h/2)), y: (Ys[i-1]+(h/2)*k1)})
        k3 = df.evalf(subs={x: (Xs[i-1]+(h/2)), y: (Ys[i-1]+(h/2)*k2)})
        k4 = df.evalf(subs={x: (Xs[i-1]+h), y: (Ys[i-1]+h*k3)})
        Ys.append(Ys[i-1]+(h/6)*(k1+2*k2+2*k3+k4))
        # print(Xs[i], Ys[i])
    return Xs, Ys


def runge_kutta_3(df, x0, y0, a, b, h):
    Xs = np.arange(a, b+h, h)
    Ys = [y0]
    # print(Xs[0], Ys[0])
    for i in range(1, len(Xs)):
        k1 = df.evalf(subs={x: Xs[i-1], y: Ys[i-1]})
        k2 = df.evalf(subs={x: (Xs[i-1]+(h/2)), y: (Ys[i-1]+(h/2)*k1)})
        k3 = df.evalf(subs={x: (Xs[i-1]+h), y: (Ys[i-1]-h*k1+2*h*k2)})
        Ys.append(Ys[i-1]+(h/6)*(k1+4*k2+k3))
        # print(Xs[i], Ys[i])
    return Xs, Ys


def runge_kutta(df, x0, y0, a, b, h):
    Xs3, Ys3 = runge_kutta_3(df, x0, y0, a, b, h)
    Xs4, Ys4 = runge_kutta_4(df, x0, y0, a, b, h)

    Rs = list(map(lambda r4, r3: (r4-r3)/(2**(4)-1), Ys4, Ys3))

    fig, ax = plt.subplots()

    ax.plot(Xs3, Ys3, 'o-')
    ax.plot(Xs4, Ys4, 's-')

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.show()
    print("-"*50)
    print("Метод Рунге-Кутта третього порядку\n")
    [print(k) for k in list(zip(Xs3, Ys3))]
    print("-"*50)
    print("Метод Рунге-Кутта четвертого порядку\n")
    [print(k) for k in list(zip(Xs4, Ys4))]
    print("-"*50)
    print("Похибка за правилом Рунге\n")
    [print(k) for k in list(zip(Xs4, Rs))]


def main():
    df = (1/(1+(x**3)*y))+2*y
    x0 = 1.5
    y0 = 2.1
    a = 1.5
    b = 2
    h = 0.05
    runge_kutta(df, x0, y0, a, b, h)


if __name__ == "__main__":
    main()
