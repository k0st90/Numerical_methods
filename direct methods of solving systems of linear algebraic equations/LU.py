import pprint
import numpy as np
from scipy import linalg

A = np.array([[0.64, 0.72,-0.83, 4.2],
            [0.58, -0.83, 1.43, -0.62],
            [0.86, 0.77, -1.83, 0.88],
            [1.32, -0.52, -0.65, 1.22]])#Матриця завдання 3.1, 3.2

B=np.array([2.23, 1.71, -0.54, 0.65])

M_33 = np.array([[1.1*np.cos(np.pi/9),np.sin(np.pi/9),-1],
                [1.2*np.cos(np.pi/4),np.sin(np.pi/4),-1],
                [2.02*np.cos(np.pi/3),np.sin(np.pi/3),-1]])#Матриця завдання 3.4

S_33 = np.array([np.power(1.1,2),np.power(1.2,2),np.power(2.02,2)])

def lu(A)->(np.ndarray,np.ndarray):
    n=len(A[0])

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

    return (L,U)



print("-"*40)
print("A:")
print(A)
L, U = lu(A)#LU розклад
print("-"*40)
print("L:")
print(L)
print("-"*40)
print("U:")
print(U)
print("-"*40)
L_1=linalg.inv(L)#Обернена матриця L
U_1=linalg.inv(U)#Обернена матриця U
A_1=np.matmul(U_1,L_1)#Обернена матриця A через LU розклад
Y=np.matmul(L_1,B)#розвʼязки рівння Ly=B
"""print("L_1:")
print(L_1)
print("-"*40)
print("U_1:")
print(U_1)
print("-"*40)
print("Y:")
print(Y)
print("-"*40)"""
X=np.matmul(U_1,Y)#розвʼязки рівння Ux=Y
print("X:")
print(X)
print("-"*40)

print("A^-1:")
print(A_1)

print("-"*40)
f=lambda x: 1 if x==-1 else U[x][x]*(f(x-1))#Детермінант
det=f(len(A[0])-1)

print("Детермінант: ",det)
print("-"*40)


print("\tЗАВДАННЯ 3.4")
KL, KU = lu(M_33)#LU розклад
KL_1=linalg.inv(KL)#Обернена матриця L
KU_1=linalg.inv(KU)#Обернена матриця U
KY=np.matmul(KL_1,S_33)#розвʼязки рівння Ly=B
K=np.matmul(KU_1,KY)#розвʼязки рівння Ux=Y
print("K=",K)
a1=K[0]/2
a3=K[1]/2*a1
a2=np.sqrt(np.power(a3,2)+np.power(a3,2)-K[2])
print("a1 =",a1,"  a2 =",a2,"  a3 =",a3)
