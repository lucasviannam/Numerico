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

def main():
    M = np.array([[1.0,1.0,1.0,3.0],[2.0,4.0,8.0,14.0],[5.0,-5.0,9.0,9.0]])
    A = np.array([[1.0,1.0,1.0],[2.0,4.0,8.0],[5.0,-5.0,9.0]])
    B = np.array([3.0,14.0,9.0])
    A1 = np.array([[1.0,1.0,1.0],[2.0,4.0,8.0],[5.0,-5.0,9.0]])
    B1 = np.array([3.0,14.0,9.0])
    Solucao1 = EliminacaoDeGauss(A,B)
    Solucao2 = np.linalg.solve(A1,B1)
    Solucao3 = EliminacaoGauss(M)
    print(Solucao1)
    print(Solucao2)
    print(Solucao3)

main()