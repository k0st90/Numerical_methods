import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x, y = sp.symbols('x y')


def euler(df, x0, y0, a, b, h):
    Xs = np.arange(a, b+h, h)
    Ys = [y0]
    print(Xs[0], Ys[0])
    for i in range(1, len(Xs)):
        Ys.append(Ys[i-1]+h*df.evalf(subs={x: Xs[i-1], y: Ys[i-1]}))
        print(Xs[i], Ys[i])

    fig, ax = plt.subplots()

    ax.plot(Xs, Ys, 'o')
    ax.plot(Xs, Ys, '-')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #plt.show()


def main():
    df = ((sp.exp(2*x)*(2+3*sp.cos(x)))/(2*y))-((3*y*sp.cos(x))/2)
    x0 = 0.2
    y0 = 0.5
    a = 0.2
    b = 1.2
    h = 0.1
    euler(df, x0, y0, a, b, h)


if __name__ == "__main__":
    print("Метод Ейлера\n")
    main()
