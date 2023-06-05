import numpy as np
"""----------------------------------------------"""
f = lambda x: (((np.log(x))**2)-(np.log(x))-2)
g = lambda x: (((np.log(x))**2)+(2*np.log(x))+1)
"""----------------------------------------------"""


def my_bisection(f, a, b, tol , k=0.3):
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("Скаляри a та b можуть не обмежувати корінь, або можуть обмежувати більше ніж один корінь!")
    m = (a + b)/2
    if np.absolute(b-a) <= 2*tol:
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        print(m,"\t", abs(m-7.3890557417868))

        return my_bisection(f, m, b, tol,k)
    elif np.sign(f(b)) == np.sign(f(m)):
        print(m,"\t", abs(m-7.3890557417868))
        return my_bisection(f, a, m, tol,k)

def chord(fun,x0, x1, e):
    condition = True
    x2 = x0 - (fun(x0)) / (fun(x1) - fun(x0)) * (x1 - x0)
    if fun(x2) * fun(x1) < 0:
        while condition:
            x3 = x2 - (fun(x2)) / (fun(x1) - fun(x2)) * (x1 - x2)
            print(x2,"\t", abs(x2-0.367879434317))
            condition = abs(x3 - x2) > e
            x2 = x3

    if fun(x2) * fun(x0) < 0:
        while condition:
            x3 = x0 - (fun(x0)) / (fun(x2) - fun(x0)) * (x2 - x0)
            print(x2,"\t", abs(x2-0.367879434317))
            condition = abs(x3 - x2) > e
            x2 = x3

    return x2


def secant(fun,x0,x1,e):
    condition = True
    while condition:
        if fun(x0) == fun(x1):
            print('Ділення на нуль!')
            break
        x2 = x0 - (x1-x0)*fun(x0)/( fun(x1) - fun(x0) )
        x0 = x1
        x1 = x2
        print(x1,"\t", abs(x1-0.367879434317))
        condition = np.absolute(x1-x0) > e
    return x2





print("----БІСЕКЦІЇ----")
"""r04 = my_bisection(f, 0.25, 1, 1e-4)
print("Корінь1 функції f за методом БІСЕКЦІЇ з точністю 10^-4:",r04)"""
"""r05 = my_bisection(f, 0.25, 1, 1e-5)
print("Корінь1 функції f за методом БІСЕКЦІЇ з точністю 10^-5: ", r05)"""
"""r06 = my_bisection(f, 0.25, 1, 1e-6)
print("Корінь1 функції f за методом БІСЕКЦІЇ з точністю 10^-6: ", r06)"""
print("-"*20)
r24 = my_bisection(f, 5, 100, 1e-4)
print("Корінь2 функції f за методом БІСЕКЦІЇ з точністю 10^-4:",r24)
"""r25 = my_bisection(f, 6, 8, 1e-5)
print("Корінь2 функції f за методом БІСЕКЦІЇ з точністю 10^-5: ", r25)
r26 = my_bisection(f, 6, 8, 1e-6)
print("Корінь2 функції f за методом БІСЕКЦІЇ з точністю 10^-6: ", r26)
print("-"*70)"""
"""print("\n\n----ХОРД----")
ch04=chord(f,0.25, 1,1e-4)
print("Корінь1 функції f за методом ХОРД з точністю 10^-4:",ch04)"""
"""ch05=chord(f,0.25, 1,1e-5)
print("Корінь1 функції f за методом ХОРД з точністю 10^-5:",ch05)
ch06=chord(f,0.25, 1,1e-6)
print("Корінь1 функції f за методом ХОРД з точністю 10^-6:",ch06)
print("-"*20)
ch24=chord(f,6,8,1e-4)
print("Корінь2 функції f за методом ХОРД з точністю 10^-4:",ch24)
ch25=chord(f,6,8,1e-5)
print("Корінь2 функції f за методом ХОРД з точністю 10^-5:",ch25)
ch26=chord(f,6,8,1e-6)
print("Корінь2 функції f за методом ХОРД з точністю 10^-6:",ch26)
print("-"*70)
"""
#print("\n\n----СІЧНИХ----")
"""se04=secant(f,0.25,0.5,1e-4)
print("Корінь1 функції f за методом СІЧНИХ з точністю 10^-4:",se04)"""
"""se05=secant(f,0.25,0.5,1e-5)
print("Корінь1 функції f за методом СІЧНИХ з точністю 10^-5:",se05)
se06=secant(f,0.25,0.5,1e-6)
print("Корінь1 функції f за методом СІЧНИХ з точністю 10^-6:",se06)
print("-"*20)
se24=secant(f,float(6),float(8),1e-4)
print("Корінь2 функції f за методом СІЧНИХ з точністю 10^-4:",se24)
se25=secant(f,float(6),float(8),1e-5)
print("Корінь2 функції f за методом СІЧНИХ з точністю 10^-5:",se25)
se26=secant(f,float(6),float(8),1e-6)
print("Корінь2 функції f за методом СІЧНИХ з точністю 10^-6:",se26)
print("-"*20)
se14=secant(g,0.25,0.5,1e-4)
print("Корінь функції g за методом СІЧНИХ з точністю 10^-4:",se14)
se15=secant(g,0.25,0.5,1e-6)
print("Корінь функції g за методом СІЧНИХ з точністю 10^-5:",se15)
se16=secant(g,0.25,0.5,1e-7)
print("Корінь функції g за методом СІЧНИХ з точністю 10^-6:",se16)"""
