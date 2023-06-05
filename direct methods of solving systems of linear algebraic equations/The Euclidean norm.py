import numpy as np
import pprint
from scipy import linalg
import matplotlib.pyplot as plt

f = lambda i,j: (1/(np.sqrt(((0.8*i*j)**2)+0.58*0.8*i*j)))
g = lambda x: x**2


def eu_norm(M)->float:
    norm=0
    n=len(M[0])
    for i in range(n):
        norm+=sum(list(map(g,M[i])))
    return np.sqrt(norm)

def cond(M)->float:
    M_1=linalg.inv(M)
    return np.multiply(eu_norm(M),eu_norm(M_1))

def A_gen(n):
    M: np.ndarray = np.array([[0.0] * n for i in range(n)])
    for i in range(n):
        for j in range(n):
            M[i][j]=f(i+1,j+1)
    return M
print("A: ")
A=A_gen(6)
print(A)

print("Евклідова норма: ",eu_norm(A))

print("число обумовленості A: ")
print(cond(A))

ns = list(range(3,9))
norms=[eu_norm(i) for i in [A_gen(n) for n in range(3,9)]]
plt.plot(ns, norms)
for i, j in zip(ns, norms):
    plt.text(i, j, str(j))
plt.xlabel('n')
plt.ylabel('Norm')
plt.title('Графік залежності норми')
plt.show()

conds=[cond(i) for i in [A_gen(n) for n in range(3,9)]]
plt.plot(ns, conds)
for i, j in zip(ns, conds):
    plt.text(i, j, str(j))
plt.xlabel('n')
plt.ylabel('cond')
plt.yscale('log')
plt.title('Графік залежності числа обумовленості')
plt.show()
