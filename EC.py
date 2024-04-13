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




def f(t,y,m): #funcao de discretizacao
    #y[0-2] Rx     y[3-5] Ry      y[6-8] Vx       y[9-11] Vy
    f = np.zeros(12,dtype=float)
    rab = distanciaModulo(y[0],y[3],y[1],y[4])
    rac = distanciaModulo(y[0],y[3],y[2],y[5])
    rbc = distanciaModulo(y[1],y[4],y[2],y[5])
    G = 1
    f[6] = G*m*((y[1]-y[0])/rab**3+(y[2]-y[0])/rac**3)
    f[7] = G*m*((y[0]-y[1])/rab**3+(y[2]-y[1])/rbc**3)
    f[8] = G*m*((y[0]-y[2])/rac**3+(y[1]-y[2])/rbc**3)
    f[9] = G*m*((y[4]-y[3])/rab**3+(y[5]-y[3])/rac**3)
    f[10]= G*m*((y[3]-y[4])/rab**3+(y[5]-y[4])/rbc**3)
    f[11]= G*m*((y[3]-y[5])/rac**3+(y[4]-y[5])/rbc**3)
    f[0] = y[6]
    f[1] = y[7]
    f[2] = y[8]
    f[3] = y[9]
    f[4] = y[10]
    f[5] = y[11]
    return f

def RK(t,yk,h,m):
    k1 = f(t,yk,m)
    k2 = f(t+h/2,yk+k1*h/2,m)
    k3 = f(t+h/2,yk+k2*h/2,m)
    k4 = f(t+h,yk+k3*h,m)
    ykplus1 = yk + h*(k1+2*k2+2*k3+k4)/6
    return ykplus1


def energia(Va,Vb,Vc,Rab,Rac,Rbc,m):
    E = m/2*(Va**2+Vb**2+Vc**2)-m**2*(1/Rab+1/Rac+1/Rbc)
    return E

def modulo(x,y): #check
    mod = math.sqrt((x**2+y**2))
    return mod

def distanciaModulo(x1,y1,x2,y2): #check
    distMod = math.sqrt((x1-x2)**2+(y1-y2)**2)
    return distMod

def prodVetorial(Rx,Ry,Vx,Vy,a,i): #check
    prod = Rx[a][i]*Vy[a][i] - Vx[a][i]*Ry[a][i]
    return prod


def momentoAngular(Rx,Ry,Vx,Vy,i,m): #check
    L = m*(prodVetorial(Rx,Ry,Vx,Vy,0,i)+prodVetorial(Rx,Ry,Vx,Vy,1,i)+prodVetorial(Rx,Ry,Vx,Vy,2,i))
    return L

def momentoLinear(Va,Vb,Vc,m): #check
    P = (Va+Vb+Vc)*m
    return P

def main():


    #criando de variáveis

    tk  = 0.0           #tempo discreto
    t0  = 0.0           #tempo inicial
    tf  = 6.325           #tempo final
    n   = 40           #numero de passos
    h   = tf/n           #tamanho do passo
    m   = 1             #massa dos corpos
    #Vx0 = 0.3471128135672417*2
    #Vy0 = 0.532726851767674*2
    Vx0 = 0.3471128135672417*2
    Vy0 = 0.532726851767674*2
    Rx  = np.zeros((3,n),dtype=float)   #primeiro indice indica qual corpo, segundo indica o k
    Ry  = np.zeros((3,n),dtype=float)
    Vx  = np.zeros((3,n),dtype=float)
    Vy  = np.zeros((3,n),dtype=float)
    E   = np.zeros(n, dtype=float)      #indice indica o k
    L   = np.zeros(n, dtype=float)
    Px  = np.zeros(n, dtype=float)
    Py  = np.zeros(n, dtype=float)
    Vmod= np.zeros(3, float)            #indice indica qual corpo
    Rab = np.single(0)
    Rac = np.single(0)
    Rbc = np.single(0)
    Yk  = np.zeros(12,dtype=float)      #variavel que recebe o resultado do método RK no passo k
    Tempos = np.linspace(t0,tf,n)

    #condicoes iniciais
    #corpo A
    Rx[0][0]  = np.single(0)
    Ry[0][0]  = np.single(0)
    Vx[0][0]  = np.single(-Vx0)
    Vy[0][0]  = np.single(-Vy0)
    #corpo B
    Rx[1][0]  = np.single(1)
    Ry[1][0]  = np.single(0)
    Vx[1][0]  = np.single(Vx0/2)
    Vy[1][0]  = np.single(Vy0/2)
    #corpo C
    Rx[2][0]  = np.single(-1)
    Ry[2][0]  = np.single(0)
    Vx[2][0]  = np.single(Vx0/2)
    Vy[2][0]  = np.single(Vy0/2)
    #carrega Yk de k=0
    for i in range(0,3):
        Yk[i]   = Rx[i][0]
        Yk[i+3] = Ry[i][0]
        Yk[i+6] = Vx[i][0]
        Yk[i+9] = Vy[i][0]
    
    Rab = distanciaModulo(Yk[0],Yk[3],Yk[1],Yk[4]) 
    Rac = distanciaModulo(Yk[0],Yk[3],Yk[2],Yk[5])
    Rbc = distanciaModulo(Yk[1],Yk[4],Yk[2],Yk[5])
    for j in range(len(Vmod)):
        Vmod[j] = modulo(Vx[j][0],Vy[j][0])

    E[0]  = energia(Vmod[0],Vmod[1],Vmod[2],Rab,Rac,Rbc,m)
    L[0]  = momentoAngular(Rx,Ry,Vx,Vy,0,m)
    Px[0] = momentoLinear(Vx[0][0],Vx[1][0],Vx[2][0],m)
    Py[0] = momentoLinear(Vy[0][0],Vy[1][0],Vy[2][0],m)


    #Integração de passo único utilizando Runge-Kutta de ordem quatro com passo h
    for i in range(1,len(Rx[0])):
        Yk = RK(tk,Yk,h,m) #chama o médoto Kunge-Kutta
        for j in range(0,3): #salva os resultados do passo k
            Rx[j][i] = Yk[j]
            Ry[j][i] = Yk[j+3]
            Vx[j][i] = Yk[j+6]
            Vy[j][i] = Yk[j+9]
        
        #Cálculo da distância entre os corpos e do módulo da posicao
        #Necessarios para o calculo da Energia e do Momento Angular do sistema
        Rab = distanciaModulo(Yk[0],Yk[3],Yk[1],Yk[4]) 
        Rac = distanciaModulo(Yk[0],Yk[3],Yk[2],Yk[5])
        Rbc = distanciaModulo(Yk[1],Yk[4],Yk[2],Yk[5])
        for j in range(len(Vmod)):
            Vmod[j] = modulo(Vx[j][i],Vy[j][i])
        #Calculo dos parâmetros Energia, Momento Agular e Momento linear do sistema
        E[i]  = energia(Vmod[0],Vmod[1],Vmod[2],Rab,Rac,Rbc,m)
        L[i]  = momentoAngular(Rx,Ry,Vx,Vy,i,m)
        Px[i] = momentoLinear(Vx[0][i],Vx[1][i],Vx[2][i],m)
        Py[i] = momentoLinear(Vy[0][i],Vy[1][i],Vy[2][i],m)



    ##axs[0,0].legend()
    #plt.close('all')
    #fig, axs = plt.subplots(3,1)
    #axs[1].plot(Tempos,E, label="Energia")
    #axs[1].legend()
    #axs[0].plot(Tempos,L, label="Momento Angular")
    #axs[0].legend()
    #axs[2].plot(Tempos,Px, label="Momento linear X")
    #axs[2].plot(Tempos,Py, label="Momento linear Y")  
    #axs[2].legend()
    #plt.savefig("./figures/ParametrosSistema.png")
    #fig1 = plt.figure()
    #plt.plot(Rx[0],Ry[0],'*' ,label="A")
    #plt.plot(Rx[1],Ry[1],'o' ,label="B")
    #plt.plot(Rx[2],Ry[2],'P' ,label="C")
    #plt.legend()
    #plt.savefig("./figures/Lamniscata.png")
    #plt.show()
    ##################################################33
    nSpline = 5
    a,b,c,d = natural_cubic_spline(Tempos,Rx[0])   
    yteste = np.zeros((len(Rx[0])-1,nSpline))
    offset = np.linspace(0,h,nSpline)
    xteste = np.zeros((len(Rx[0])-1,nSpline))
    
    for i in range(len(Rx[0])-1):
            for j in range(nSpline):
                xteste[i][j] = Tempos[i]+offset[j]

    for i in range(len(Rx[0])-1):
        for j in range(nSpline):
            yteste[i][j] = a[i] + b[i]*(xteste[i][j]-xteste[i][0]) + c[i]*(xteste[i][j]-xteste[i][0])**2 + d[i]*(xteste[i][j]-xteste[i][0])**3

    linha = [':','--']
    for i in range(len(yteste)-1):
        plt.plot(xteste[i],yteste[i],linha[i%2])
    plt.plot(Tempos,Rx[0],'o')
    plt.xlim(2,4.5)
    plt.savefig("./figures/SplineEC.png")
    plt.show()
    
main()