import sympy as sp
x = sp.symbols('x')
f = sp.Pow(sp.sin(8 + sp.simplify(sp.log(sp.Abs(x), 10)) + sp.Pow(sp.log(sp.Abs(x), 10), 2)), 3)
X = []
Y = []
x0 = 10.9
for i in range(6):
    X.append(x0+i*0.1)
for xi in X:
    Y.append(f.evalf(subs={x: xi}))


def newton_pol():
    dif = [-0.256883680464095, -0.115384468, -0.006983865, -0.001720625, -5.1834e-06, 1.39388e-05]
    N = [f'({dif[0]})']
    for i in range(1, len(dif)):
        tempp = []
        for j in range(i):
            tempp.append(f'(x-{X[j]})')
        N.append(f"({dif[i]})*"+"*".join(tempp))
    return sp.sympify("+".join(N)).expand()


def main():
    print("X: ")
    print(X)
    print("Y: ")
    print(Y)
    N = newton_pol()
    print(N)
    print()
    print("\nNm(x1+x2): ", N.evalf(subs={x: X[1]+X[2]}))
    sp.plotting.plot(N, (x, -2, 35), markers=[{'args': [X, Y, 'ro']}])  # type: ignore


if __name__ == "__main__":
    main()
