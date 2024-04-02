import matplotlib.pyplot as plt
import math
import numpy as np

def teste(matriz):
    print(matriz)
    return

def modulo(x,y):
    mod = math.sqrt((x**2+y**2))
    return mod

def prodVetorial(Rx,Ry,Vx,Vy,a,i):
    prod = Rx[a][i]*Vy[a][i] - Vx[a][i]*Ry[a][i]
    return prod

def momentoAngular(Rx,Ry,Vx,Vy,i,m):
    L = m*(prodVetorial(Rx,Ry,Vx,Vy,0,i)+prodVetorial(Rx,Ry,Vx,Vy,1,i)+prodVetorial(Rx,Ry,Vx,Vy,2,i))
    return L

def distanciaModulo(x1,y1,x2,y2):
    distMod = math.sqrt((x1-x2)**2+(y1-y2)**2)
    return distMod

def main():

    m = 10
    t0  = 0
    tf  = 10
    n   = 500
    h   = (10+0)/500

    Rx = np.array([[-1,-0.75,-0.5],[0,0.2,0.4],[1,0.75,0.5]])
    Ry = np.array([[0,0.2,0.5],[0,-0.2,-0.4],[0,0.3,0.4]])
    Vx = np.array([[1,2,3],[2,2,2],[-1.2,-2.5,-3]])
    Vy = np.array([[1,1.2,1.4],[-1.2,-1.3,-1.4],[1.4,1.5,1.6]])
    print(Rx)
    print(Ry/2)
    teste(Rx+Ry/2)

    plt.plot(Rx[0],Ry[0])

main()