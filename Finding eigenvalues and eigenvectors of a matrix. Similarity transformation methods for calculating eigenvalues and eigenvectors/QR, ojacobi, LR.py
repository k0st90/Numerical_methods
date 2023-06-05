import numpy as np
import cmath


def lu(A) -> (np.ndarray, np.ndarray):
    n = len(A[0])

    L: np.ndarray = np.array([[0.0] * n for i in range(n)])
    U: np.ndarray = np.array([[0.0] * n for i in range(n)])

    for j in range(n):
        L[j][j] = 1.0

        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = A[i][j] - s1

        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (A[i][j] - s2) / U[j][j]

    return (L, U)


def householder_vector(x):
    n = len(x)
    e1 = np.zeros(n)
    e1[0] = 1
    v = x + np.sign(x[0]) * np.linalg.norm(x) * e1
    return v / np.linalg.norm(v)


def householder_matrix(v):
    n = len(v)
    return np.identity(n) - 2 * np.outer(v, v)


def QR(A):
    m, n = A.shape
    Q = np.identity(m)
    R = np.copy(A)

    for j in range(min(n, m-1)):
        v = householder_vector(R[j:, j])
        Hj = np.identity(m)
        Hj[j:, j:] = householder_matrix(v)
        Q = np.dot(Q, Hj)
        R = np.dot(Hj, R)

    return Q, R


def lr(M, eps):

    n = M.shape[0]
    RP = np.zeros((n, n))
    R = np.ones((n, n))
    M1 = np.copy(M)
    counter = 1
    while np.linalg.norm(R-RP) > eps:
        RP = np.copy(R)
        L, R = lu(M1)
        print(counter, np.diag(R))
        M1 = np.matmul(R, L)
        counter += 1
    return np.diag(M1)


def qr(M, eps):
    n = M.shape[0]
    RP = np.zeros((n, n))
    R = np.ones((n, n))
    M1 = np.copy(M)
    counter = 1
    while np.linalg.norm(R-RP) > eps:
        RP = np.copy(R)
        Q, R = QR(M1)
        print(counter, np.diag(R))
        # print(M1)
        M1 = np.matmul(R, Q)
        counter += 1
    return np.diag(M1)


def queq(S):
    a = 1
    b = -(S[0][0]+S[1][1])
    c = S[0][0]*S[1][1]-S[0][1]*S[1][0]

    d = (b**2) - (4*a*c)

    sol1 = (-b-cmath.sqrt(d))/(2*a)
    sol2 = (-b+cmath.sqrt(d))/(2*a)

    return sol1, sol2


def qr_comp(M, eps):
    n = M.shape[0]
    if n != 4:
        return
    R = np.ones((n, n))
    M1 = np.copy(M)
    counter = 1
    mask = np.zeros(M1.shape, dtype=bool)
    mask[2][0] = 1
    mask[3][0] = 1
    mask[3][1] = 1

    while M1[mask].all() > eps:
        Q, R = np.linalg.qr(M1)
        M1 = np.matmul(R, Q)
        counter += 1
    if np.abs(M1[1][0]) > eps and np.abs(M1[3][2]) < eps:
        e1, e2 = queq(M1[0:2, 0:2])
        e3 = M1[2][2]
        e4 = M1[3][3]
        return [e1, e2, e3, e4]

    if np.abs(M1[3][2]) > eps and np.abs(M1[1][0]) < eps:
        e1, e2 = queq(M1[2:4, 2:4])
        e3 = M1[0][0]
        e4 = M1[1][1]
        return [e1, e2, e3, e4]

    if np.abs(M1[3][2]) > eps and np.abs(M1[1][0]) > eps:
        e1, e2 = queq(M1[0:2, 0:2])
        e3, e4 = queq(M1[2:4, 2:4])
        return [e1, e2, e3, e4]

    if np.abs(M1[2][1]) > eps:
        e1, e2 = queq(M1[1:3, 1:3])
        e3 = M1[0][0]
        e4 = M1[3][3]
        return [e1, e2, e3, e4]


def find_max(M):
    n = M.shape[0]
    max_el = [0, 0, 0]
    for i in range(n):
        for j in range(n):
            if np.abs(M[i][j]) > np.abs(max_el[0]) and i < j:
                max_el[0] = M[i][j]
                max_el[1] = i
                max_el[2] = j

    return max_el


def rotation(M):
    n = M.shape[0]
    a, l, m = find_max(M)
    if M[l][l] == M[m][m]:
        res = np.pi/4
    else:
        res = 0.5 * np.arctan((2*a)/(M[l][l]-M[m][m]))
    H = np.eye(n)
    H[l][l] = np.cos(res)
    H[l][m] = -np.sin(res)
    H[m][l] = np.sin(res)
    H[m][m] = np.cos(res)
    return H


def ojacobi(M, eps):
    counter = 1
    A = np.copy(M)
    cond = np.inf
    while cond > eps:
        H = rotation(A)
        A = np.dot(np.dot(np.transpose(H), A), H)
        print(counter, np.diag(A))
        cond = np.sqrt(np.sum(np.square(np.tril(A, -1))))
        counter += 1
    return np.diag(A)


def main():

    A = np.array([[2., 1.5, 3.5, 4.5],
                  [1.5, 2., 2., 1.6],
                  [3.5, 2., 2., 1.7],
                  [4.5, 1.6, 1.7, 2.]])

    A63 = np.array([[0.72, 3.54, 7.28, 0.33],
                    [-0.28, -0.72, 3.04, 0.22],
                    [1., 0.35, -0.78, 1.],
                    [7.03, -5.04, -3.75, 3.41]])

    print("\n\tLR метод з точністю 0.01")
    print(lr(A, 0.01))
    print("\n\tLR метод з точністю 0.001")
    print(lr(A, 0.01))
    print("\n\tQR метод з точністю 0.01")
    print(qr(A, 0.01))
    print("\n\tQR метод з точністю 0.001")
    print(qr(A, 0.001))
    print("\n\tМетод обертань Якобі з точністю 0.01")
    print(ojacobi(A, 0.01))
    print("\n\tМетод обертань Якобі з точністю 0.001")
    print(ojacobi(A, 0.001))
    print("-"*50)
    print("\n\tQR метод завдання 6.3 з точністю 0.01")
    print(qr_comp(A63, 0.01))
    print("\n\tQR метод завдання 6.3 з точністю 0.001")
    print(qr_comp(A63, 0.001))
    print("\nВласні значення обраховані за допомогою математичного пакета:\n", *np.linalg.eig(A63)[0])


if __name__ == '__main__':
    main()
