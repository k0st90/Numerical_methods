import numpy as np
from sympy import *

"""----------------------------------------------"""
x = Symbol("x")
f = (((ln(x))**2)-(ln(x))-2)
df=diff(f,x)
g = (((ln(x))**2)+(2*ln(x))+1)
dg=diff(g,x)
p=lambda x: (-0.25*np.exp(x))+((x**2)/2)+0.5
q=lambda x: np.cbrt(0.1*(x**2)-(0.4*x)-2)
"""----------------------------------------------"""

def newton(f, df, x0, tol):
    print(x0,"\t",abs(x0-7.3890560989305873374))
    if abs(f.evalf(30,subs={x:x0})/df.evalf(30,subs={x:x0})) < tol:
        return x0
    else:

        return newton(f, df, x0 - f.evalf(30,subs={x:x0})/df.evalf(30,subs={x:x0}), tol)

fff=df.evalf(30,subs={x:8})
def s_newton(f, x0, tol,fff):
    print(x0,"\t",abs(x0-7.3890560989305873374))
    if abs((f.evalf(30,subs={x:x0}))/fff) < tol:
        return x0
    else:
        #print(x0)
        return s_newton(f, x0 - f.evalf(30,subs={x:x0})/fff, tol,fff)

def simple_iter(f,a,b,tol):
    counter=0;
    x=(a+b)/2
    print(x)
    z=x
    x=f(x)
    while ((abs(x-z))>=tol):
        print(x,abs(x-0.2133086597637499))
        z=x
        x=f(x)
        #print(counter,x)
        counter+=1
    return x

"""print("----НЬЮТОНА----")
r04 = newton(f, df, 0.25, 0.0001)
print("Корінь1 функції f за методом НЬЮТОНА з точністю 10^-4:",r04)
r14 = newton(f, df, 0.25, 0.00001)
print("Корінь1 функції f за методом НЬЮТОНА з точністю 10^-5:",r14)
r24 = newton(f, df, 0.25, 0.000001)
print("Корінь1 функції f за методом НЬЮТОНА з точністю 10^-6:",r24)
print("-"*20)
r024 = newton(f, df, 8, 0.0001)
print("Корінь2 функції f за методом НЬЮТОНА з точністю 10^-4:",r024)
r124 = newton(f, df, 8, 0.00001)
print("Корінь2 функції f за методом НЬЮТОНА з точністю 10^-5:",r124)
r224 = newton(f, df, 8, 0.000001)
print("Корінь2 функції f за методом НЬЮТОНА з точністю 10^-6:",r224)
print("-"*20)
print("----CПРОЩЕНИЙ МЕТОД НЬЮТОНА----")
r230 = s_newton(f, 0.1, 0.0001)
print("Корінь1 функції f за методом спрощеного ньютона з точністю 10^-4:",r230)
r231 = s_newton(f, 0.1, 0.00001)
print("Корінь1 функції f за методом спрощеного ньютона з точністю 10^-5:",r231)
r232 = s_newton(f, 0.1, 0.000001)
print("Корінь1 функції f за методом спрощеного ньютона з точністю 10^-6:",r232)
print("-"*20)
r240 = s_newton(f, 8, 0.0001,fff)
print("Корінь2 функції f за методом спрощеного ньютона з точністю 10^-4:",r240)
r241 = s_newton(f, 8, 0.00001)
print("Корінь2 функції f за методом спрощеного ньютона з точністю 10^-5:",r241)
r242 = s_newton(f, 8, 0.000001)
print("Корінь2 функції f за методом спрощеного ньютона з точністю 10^-6:",r242)
print("-"*20)"""

print("----ПРОСТОЇ ІТЕРАЦІЇ----")
"""r250 = simple_iter(p,0, 0.5, 0.0001)
print("Корінь функції p за методом простої ітерації з точністю 10^-4:",r250)
r251 = simple_iter(p,0, 0.5, 0.00001)
print("Корінь функції p за методом простої ітерації з точністю 10^-5:",r251)"""
r252 = simple_iter(p,0, 0.5, 0.000001)
print("Корінь функції p за методом простої ітерації з точністю 10^-6:",r252)
print("-"*20)
"""r350 = simple_iter(q,-1.2, -1, 0.0001)
print("Корінь функції q за методом простої ітерації з точністю 10^-4:",r350)
r351 = simple_iter(q,-1.2, -1, 0.00001)
print("Корінь функції q за методом простої ітерації з точністю 10^-5:",r351)
r352 = simple_iter(q,-1.2, -1, 0.000001)
print("Корінь функції q за методом простої ітерації з точністю 10^-6:",r352)
print("-"*20)"""
