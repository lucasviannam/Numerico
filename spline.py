import matplotlib.pyplot as plt
import math
import numpy as np

def EliminacaoGauss(M):
    n = len(M)
    Solucao = np.zeros(n)
    for i in range(n-1):
        for j in range(i+1,n):
            M[j][i] = M[j][i]/M[i][i]
            for k in range(i+1,n+1):
                M[j][k] = M[j][k] - M[j][i]*M[i][k]
    for i in range(n-1,-1,-1):
        for j in range(i+1,n):
            M[i][n]=M[i][n]-M[i][j]
        M[i][i] = M[i][n]/M[i][i]
        for j in range(i-1,-1,-1):
            M[j][i] = M[j][i]*M[i][i]
    for i in range(n):
        Solucao[i] = M[i][i]
    return Solucao



def natural_cubic_spline(x, y):
    n = len(x) - 1
    h = x[1] - x[0]

    A = np.zeros((n,n),dtype=float)
    A[0][0] = 1.0
    A[n-1][n-1] = 1.0
    for i in range(1,n-1):
        A[i][i-1] = h
        A[i][i] = 4*h
        A[i][i+1] = h
    B = np.zeros(n)
    for i in range(1,n-1):
        B[i] = (y[i+1]-2*y[i]+y[i-1])*3/h
    c = np.zeros(n)
    c = np.linalg.solve(A,B)
    b = np.zeros(n)
    d = np.zeros(n)
    for i in range(n-1):
        b[i] = (y[i+1]-y[i])/h - h*(c[i+1]+2*c[i])/3
        d[i] = (c[i+1]-c[i])/(3*h)

    return(y,b,c,d)




def main():
    #M = np.array([[1.0,1.0,1.0,3.0],[2.0,4.0,8.0,14.0],[5.0,-5.0,9.0,9.0]])
    #print(len(M))
    #print(M[1][0])
    #Solucao = EliminacaoGauss(M)
    Npontos = 20
    y = np.zeros(Npontos)
    x = np.linspace(1,20,Npontos)

    for i in range(Npontos):
        y[i] = 2*x[i]**3-5*x[i]**2+3*x[i]+15
    #print(x)
    #print(y)
    a,b,c,d = natural_cubic_spline(x,y)
    yteste = np.zeros((Npontos-1,20))
    offset = np.linspace(0,1,20)
    xteste = np.zeros((Npontos-1,20))
    for i in range(Npontos-1):
            for j in range(20):
                xteste[i][j] = x[i]+offset[j]
    for i in range(Npontos-1):
        for j in range(20):
            yteste[i][j] = a[i] + b[i]*(xteste[i][j]-xteste[i][0]) + c[i]*(xteste[i][j]-xteste[i][0])**2 + d[i]*(xteste[i][j]-xteste[i][0])**3

    print(a)
    print(b)
    print(c)
    print(d)
    plt.plot(x,y,'o')
    for i in range(len(yteste)):
        plt.plot(xteste[i],yteste[i])
    #plt.ylim(0,8000)
    plt.show()
    #coefficients = natural_cubic_spline(x, y)
    #for i, coef in enumerate(coefficients):
    #    print(f"S{i}(x) = {coef[0]} + {coef[1]}(x - {x[i]}) + {coef[2]}(x - {x[i]})^2 + {coef[3]}(x - {x[i]})^3")
    #plt.plot(x,y,'o')
    #plt.show()

    #print(Solucao)

main()

