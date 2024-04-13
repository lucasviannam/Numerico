import matplotlib.pyplot as plt
import math
import numpy as np

def EliminacaoDeGauss(A,B):
    n =  len(A)
    mult = np.zeros(n-1)
    Solucao = np.zeros(n)
    for i in range(n-1):
        for j in range(i+1,n):
            mult[j-1] = A[j][i]/A[i][i]
            A[j][i+1:] = A[j][i+1:] - A[i][i+1:]*mult[j-1]
            B[j] = B[j] - B[i]*mult[j-1]
            A[j][i] = 0.0
    for i in range(n-1,-1,-1):
        B[i] = B[i]-np.sum(A[i][i+1:])
        A[i][i] = B[i]/A[i][i]
        A[:i,i] = A[:i,i]*A[i][i]
    for i in range(n):
        Solucao[i] = A[i][i]
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
    #c = np.linalg.solve(A,B)
    c = EliminacaoDeGauss(A,B)
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
    yexato = np.zeros(45)
    xexato = np.linspace(0,7,45)
    yteste = np.zeros((Npontos-1,20))
    offset = np.linspace(0,1,20)
    xteste = np.zeros((Npontos-1,20))
    for i in range(Npontos-1):
            for j in range(20):
                xteste[i][j] = x[i]+offset[j]
    for i in range(Npontos-1):
        for j in range(20):
            yteste[i][j] = a[i] + b[i]*(xteste[i][j]-xteste[i][0]) + c[i]*(xteste[i][j]-xteste[i][0])**2 + d[i]*(xteste[i][j]-xteste[i][0])**3
    for i in range(len(yexato)):
        yexato[i] = 2*xexato[i]**3-5*xexato[i]**2+3*xexato[i]+15 


    print(a)
    print(b)
    print(c)
    print(d)
    plt.plot(x[:Npontos-1],y[:Npontos-1],'o',label="Pontos amostrados")
    linha = [':','--']
    for i in range(len(yteste)-1):
        plt.plot(xteste[i],yteste[i],linha[i%2])
    plt.legend()
    plt.savefig("./figures/SplineTeste.png")
    plt.plot(xexato,yexato,'k+',label="Pontos não amostrados")
    plt.xlim(1,7)
    plt.ylim(0,500)
    plt.legend()
    plt.savefig("./figures/SplineTesteDetalhe.png")
    plt.show()

main()

