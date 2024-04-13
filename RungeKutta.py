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
    k1 = f(t,yk)
    k2 = f(t+h/2,yk+k1*h/2)
    k3 = f(t+h/2,yk+k2*h/2)
    k4 = f(t+h,yk+k3*h)
    ykplus1 = yk + h*(k1+2*k2+2*k3+k4)/6
    return ykplus1

def Yteorico(C,V,lambda1,lambda2,lambda3,t):
    Y = C[0]*np.exp(lambda1*t)*V[:,0] + C[1]*np.exp(lambda2*t)*V[:,1] + C[2]*np.exp(lambda3*t)*V[:,2]
    return Y
def norma(A,B):
    if len(A) != len(B):
        print("vetores de tamanhos diferentes")
        return
    soma = 0.0
    for i in range(len(A)):
        soma = soma + (A[i]-B[i])**2
    soma = np.sqrt(soma)
    return soma

def main():
    #criando de variáveis
    t0  = 0.0           #tempo inicial
    tf  = 1.0           #tempo final
    m = 8               #numero de casos
    nmax = 16*2**(m-1)
    tk  = np.zeros(nmax)
    h   = np.zeros(m,dtype=float)


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
    Y = np.zeros([3,nmax],dtype=float)
    Y[0][0] = 3
    Y[1][0] = 4
    Y[2][0] = 6
    C = np.zeros(3)
    C = np.matmul(Vlinha,Y[:,0])
    Yexato = np.zeros((3,nmax))
    #Integração de passo único utilizando Runge-Kutta de ordem quatro com passo h
    Yk = np.zeros((3,nmax))
    YfinalAprox = np.zeros((3,m))
    Yk = Y[:,0]
    for j in range(1,m+1):
        n = 16*2**(j-1)
        tk, h[j-1]= np.linspace(t0,tf,n,retstep=True)

        for i in range(1,n):
            Yk = RK(tk[i],Y[:,i-1],h[j-1])      #chama o médoto Kunge-Kutta
            Y[:,i] = Yk

        #Tabela de convergencia
        YfinalAprox[:,j-1] = Yk
        e=p=q=r=0;
        if j>1:
            Yfinal = Yteorico(C,V,lambda1,lambda2,lambda3,tk[n-1])
            q = abs(norma(Yfinal,YfinalAprox[:,j-2])/norma(Yfinal,YfinalAprox[:,j-1]));
            r = h[j-2]/h[j-1];
            p = math.log(q)/math.log(r);
            e = norma(Yfinal,YfinalAprox[:,j-2])
        print("%5d & %9.3e & %9.3e & %9.3e \\\\" % (n,h[j-1],e,p));

        #Plot
        string = "n = " + str(n)
        linha = ['--','-.',':','.',',','-']
        if (j+1)%2==0: 
            plt.plot(tk[:n-1],Y[0][:n-1],linha[np.int16(j/2)],label=string,)
    
    for i in range(nmax):
        #Yexato[:,i] = C[0]*np.exp(lambda1*tk[i])*V[:,0] + C[1]*np.exp(lambda2*tk[i])*V[:,1] + C[2]*np.exp(lambda3*tk[i])*V[:,2]
        Yexato[:,i] = Yteorico(C,V,lambda1,lambda2,lambda3,tk[i])
    #print(Y)
    plt.plot(tk,Yexato[0],label="Teorico")
    #plt.plot(tk,Y[0],label="simulacao")
    #plt.title("Variável Y[0] simulada com diferentes tamanhos de passo")
    plt.legend()
    plt.savefig("./figures/RungeKuttaLonge.png")
    #plt.title("Variável Y[0] simulada com diferentes tamanhos de passo \n melhor visualização")
    plt.xlim(0.792,0.807)
    plt.ylim(-1400,-1200)
    plt.savefig("./figures/RungeKuttaPerto.png")
    plt.show()
    #fig1 = plt.figure()
    #plt.legend()
    #plt.show()
    print("fim")
main()