import numpy as np
import sympy as sp

x = sp.Symbol("x")
y = sp.Symbol("y")


def plotcircles(M):
    global x
    global y
    n = M.shape[0]
    D = np.diag(np.diag(M))
    U = np.triu(M, 1)
    L = np.tril(M, -1)
    k = np.repeat([1.], n)

    for i in range(n):
        k[i] = np.sum(U[i]) + np.sum(L[i])

    for i in range(n):
        print(f"Інтервал власного значення {i+1}: [{(np.sum(D[0])-k[i])}, {np.sum(D[0])+k[i]}]")

    # print(k)
    p1 = sp.plot_parametric(
        (np.sum(D[0])+k[0]*sp.cos(x), k[0]*sp.sin(x)), (x, -5, 5), show=False)
    for j in range(1, n):
        p2 = sp.plot_parametric(
            (np.sum(D[j])+k[j]*sp.cos(x), k[j]*sp.sin(x)), (x, -5, 5), show=False)
        p1.append(p2[0])
    p1.show()


def main():

    A = np.array([[2, 1.5, 3.5, 4.5],
                  [1.5, 2., 2., 1.6],
                  [3.5, 2., 2., 1.7],
                  [4.5, 1.6, 1.7, 2.]])
    plotcircles(A)


if __name__ == '__main__':
    main()
