import numpy as np
from colorama import Fore


def powermethod(M, tol, reverse=False):
    if reverse:
        print(Fore.YELLOW + "\n\tМетод зворотних ітерацій")
        M_1 = (np.linalg.inv(M))
        M = np.copy(M_1)
    else:
        print(Fore.MAGENTA + "\n\tСтепеневий метод")
    n = M.shape[0]
    Y1 = np.ones(n)

    Za = np.zeros(n)
    Z = np.ones(n)
    la = 1
    counter = 1
    while np.linalg.norm(np.abs(Z)-np.abs(Za)) > tol:
        Za = np.copy(Z)
        Z = np.dot(M, Y1)
        la = Z[0]/Y1[0]
        if reverse:
            print(counter, " - ", 1/la)
        else:
            print(counter, " - ", la)
        Y1 = Z/np.linalg.norm(Z, np.inf)
        counter += 1
    if reverse:
        return 1/la, Y1
    return la, Y1


def complex_power_method(M, tol):

    print(Fore.CYAN + "\n\tКомплексний степеневий метод з використанням відношення Релея")

    n = M.shape[0]
    x = np.ones(n)
    y = np.ones(n)
    yp = np.zeros(n)

    eigenvalue = 0
    eigenvector = x

    counter = 1
    while np.linalg.norm(np.abs(y)-np.abs(yp)) > tol:
        yp = np.copy(y)
        y = np.dot(M, x)
        rk = np.dot(y, x) / np.dot(x, x)
        x = y / np.linalg.norm(y, np.inf)
        print(counter, " - ", eigenvalue)
        eigenvalue = np.copy(rk)
        eigenvector = np.copy(x)
        counter += 1
    return eigenvalue, eigenvector


def main():
    A = np.array([[1, 0.5, 1.2, -1],
                  [0.5, 2, -0.5, 0],
                  [1.2, -0.5, -1, 1.4],
                  [-1, 0, 1.4, 1]])  # type: ignore

    print(Fore.GREEN + str(A))
    lam, v = powermethod(A, 0.001)
    print("\nМаксимальне за абсолютною величиною власне значення:", lam,
          "\nВідповідний власний вектор: ", v)
    lambdi = np.linalg.eig(A)[0]
    lmax = np.max(np.abs(lambdi))
    llmin = np.min(np.abs(lambdi))
    print(f"Абсолютна похибка: {np.abs(lmax-np.abs(lam))}")
    print(f"Відносна похибка: {(np.abs((lmax-np.abs(lam)))/lmax)*100} %")
    # print(np.linalg.eig(A))
    # print(np.dot(A, v))
    # print(lam*v)

    lamin, vmin = powermethod(A, 0.001, True)
    print("\nМінімальне за абсолютною величиною власне значення:", lamin,
          "\nВідповідний власний вектор: ", vmin)
    print(f"Абсолютна похибка: {np.abs(llmin-np.abs(lamin))}")
    print(f"Відносна похибка: {(np.abs((llmin-np.abs(lamin)))/llmin)*100} %")
    # print(np.dot(A, vmin))
    # print(lamin*vmin)

    cla, cv = complex_power_method(A, 0.001)
    print("\nМаксимальне за абсолютною величиною власне значення:", cla,
          "\nВідповідний власний вектор: ", cv)
    print(f"Абсолютна похибка: {np.abs(lmax-np.abs(cla))}")
    print(f"Відносна похибка: {(np.abs((lmax-np.abs(cla)))/lmax)*100} %")
    # print(np.dot(A, cv))
    # print(cla*cv)


if __name__ == "__main__":
    main()
