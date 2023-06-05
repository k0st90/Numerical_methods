import numpy as np
import sympy as sp
from scipy import linalg

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


def eu_norm(M) -> float:
    """ Returns euclidean norm of matrix """
    return np.sqrt(np.sum(np.square(M)))


def cond(M) -> float:
    M_1 = linalg.inv(M)
    return np.multiply(eu_norm(M), eu_norm(M_1))


def main():

    A_k = np.array([[1, 1.5, 2.5, 3.5],  # type: ignore
                    [1.5, 1, 2, 1.6],
                    [2.5, 2, 1, 1.7],
                    [3.5, 1.6, 1.7, 1]])

    A = np.array([[1, 0.5, 1.2, -1],  # type: ignore
                  [0.5, 2, -0.5, 0],
                  [1.2, -0.5, -1, 1.4],
                  [-1, 0, 1.4, 1]])

    print("Евклідова норма: ", eu_norm(A),
          "\nЧисло обумовленості: ", cond(A))
    plotcircles(A)


if __name__ == '__main__':
    main()
