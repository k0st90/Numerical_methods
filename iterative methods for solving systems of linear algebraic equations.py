import numpy as np

A41 = np.array([[0.13, 0.27, -0.22, -0.18],
                [-0.21, 0, -0.45, 0.18],
                [0.12, 0.13, -0.33, 0.18],
                [0.33, -0.05, 0.06, -0.28]])

B41 = np.array([1.21, -0.33, -0.48, -0.17])

A42 = np.array([[4.238, 0.329, 0.256, 0.425],
                [0.249, 2.964, 0.351, 0.127],
                [0.365, 0.217, 2.897, 0.168],
                [0.178, 0.294, 0.432, 3.701]])

B42 = np.array([0.56, 0.38, 0.778, 0.749])


def dlu(M: np.ndarray):
    D = np.diag(np.diag(M))
    U = np.triu(M, 1)
    L = np.tril(M, -1)
    return D, L, U


def eu_norm(M: np.ndarray) -> float:
    """ Returns euclidean norm of matrix """
    return np.sqrt(np.sum(np.square(M)))


def alphabeta_jacobi(M: np.ndarray, B: np.ndarray):
    D, L, U = dlu(M)
    a = np.matmul((np.linalg.inv(D)*-1), np.add(L, U))
    b = np.matmul(np.linalg.inv(D), B)
    return a, b


def jacobi(A, B, e: int):
    if np.all(np.diag(np.abs(A)) > (np.sum(np.abs(A), axis=1) - np.diag(np.abs(A)))):
        a, b = alphabeta_jacobi(A, B)
        x = simp(a, b, e, 3)
        return x
    else:
        print("Метод не збігається")


def simp(A: np.ndarray, B: np.ndarray, e: float, ftv=0) -> np.ndarray:
    if eu_norm(A) < 1:
        n = len(A[0])
        Xp = np.zeros(n)
        X = np.copy(B)
        counter = 0
        while eu_norm(X-Xp) > e:
            Xp = np.copy(X)
            temp = np.zeros(n)
            for i in range(n):
                x = np.matmul(A[i], X)+B[i]
                temp[i] = x
            X = np.copy(temp)
            counter += 1
            print(X, counter)
            if counter == 3 and ftv!=0:
                print("Локальна похибка третьої ітерації: ",
                      (eu_norm(A)/(1-eu_norm(A)))*(eu_norm(X-Xp)))
        return X
    else:
        print("Метод не збігається")


def g_sei(A: np.ndarray, B: np.ndarray, e: int):
    n = len(A[0])
    a, b = alphabeta_jacobi(A, B)
    X = np.copy(b)
    Xp = np.zeros(n)
    counter = 0
    while eu_norm(X-Xp) > e:
        Xp = np.copy(X)
        for i in range(n):
            x = np.matmul(a[i], X)+b[i]
            X[i] = x
        counter += 1
        print(X, counter)
        if counter == 3:
            print("Локальна похибка третьої ітерації: ",
                  (eu_norm(dlu(A)[2])/(1-eu_norm(a)))*(eu_norm(X-Xp)))
    return X


print("-"*54)
print("Метод простої ітерації: ")
xs = simp(A41, B41, 0.001)
print("\nРозвʼязок", xs)

print("-"*54)
print("Метод Якобі: ")
xj = jacobi(A42, B42, 0.001)
print("\nРозвʼязок", xj, "\nНевʼязка розвʼязку: ", np.matmul(A42, xj)-B42)
print("-"*54)
print("Метод Гауса-Зейделя: ")
xg = g_sei(A42, B42, 0.001)
print("\nРозвʼязок", xg, "\nНевʼязка розвʼязку: ", np.matmul(A42, xg)-B42)
print("-"*54)
