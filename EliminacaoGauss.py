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



def main():
    M = np.array([[1.0,1.0,1.0,3.0],[2.0,4.0,8.0,14.0],[5.0,-5.0,9.0,9.0]])
    A = np.array([[1.0,1.0,1.0],[2.0,4.0,8.0],[5.0,-5.0,9.0]])
    B = np.array([3.0,14.0,9.0])
    Solucao1 = EliminacaoGauss(M)
    Solucao2 = np.linalg.solve(A,B)
    print(Solucao1)
    print(Solucao2)

main()