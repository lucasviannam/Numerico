import matplotlib.pyplot as plt
import math
import numpy as np

def f(t,y):         #funcao de discretizacao
    f = np.zeros(3,dtype=float)
    V = np.array([[1,4,6],[1,8,3],[1,12,2]])     
    Lambda = np.array([[8,0,0],[0,3,0],[0,0,5]]) 
    Vlinha = np.linalg.inv(V)
    A = np.matmul(V,Lambda)
    A = np.matmul(A,Vlinha)
    f = np.matmul(A,y)
    return f

def RK(t,yk,h):
    k1 = f(t,yk)*h
    k2 = f(t+h/2,yk+k1*h/2)*h
    k3 = f(t+h/2,yk+k2*h/2)*h
    k4 = f(t+h,yk+k3*h)*h
    ykplus1 = yk + h*(k1+2*k2+2*k3+k4)/6
    return ykplus1

def main():
    #criando de variáveis
    t0  = 0.0           #tempo inicial
    tf  = 0.5          #tempo final
    tk  = 0.0           #t iterador
    n   = 100           #numero de passos
    h   = (tf-t0)/n     #tamanho do passo
    tk  = np.zeros(n)
    for i in range(n):
        tk[i] = t0+i*h
    lambda1 = 8
    lambda2 = 3
    lambda3 = 5
    #Lambda1 = 8 Lambda2 = 3 Lambda3 = 5
    #V1=[1,1,1] V2=[4,8,12] V3=[6,3,2]
    V = np.array([[1,4,6],[1,8,3],[1,12,2]])     
    Lambda = np.array([[8,0,0],[0,3,0],[0,0,5]]) 
    Vlinha = np.linalg.inv(V)
    A = np.matmul(V,Lambda)
    A = np.matmul(A,Vlinha)
    Y = np.zeros([3,n],dtype=float)
    Y[0][0] = 3
    Y[1][0] = 4
    Y[2][0] = 6
    C = np.zeros(3)
    C = np.matmul(Vlinha,Y[:,0])
    Yexato = np.zeros((3,n))
    for i in range(n):
        Yexato[:,i] = C[0]*np.exp(lambda1*tk[i])*V[:,0] + C[1]*np.exp(lambda2*tk[i])*V[:,1] + C[2]*np.exp(lambda3*tk[i])*V[:,2]
    #Integração de passo único utilizando Runge-Kutta de ordem quatro com passo h
    Yk = np.zeros((3,n))
    Yk = Y[:,0]
    for i in range(1,n):
        Yk = RK(tk[i],Y[:,i-1],h)      #chama o médoto Kunge-Kutta
        Y[:,i] = Yk
    #print(Y)
    #plt.plot(tk,Yexato[0],label="solucao")
    plt.plot(tk,Y[0],label="simulacao")
    plt.legend()
    plt.show()
    #fig1 = plt.figure()
    #plt.legend()
    #plt.show()
    print("fim")
main()