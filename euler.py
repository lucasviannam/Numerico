import matplotlib.pyplot as plt
import math
import numpy as np

def f(t):
    ''' Definicao da EDO'''
    f = math.exp(t)*(2*math.cos(2*t) + math.sin(2*t))
    return f

#-----------------------------------------------------------------------------

def euler(t,tf,n):
    '''Metodo de Euler'''
    y = [0]*(n+1)
    tempos = [0]*(n+1)
    y[0] = 1
    tempos[0] = t
    delta_t = (tf-t)/n
    for i in range(0,n):
        y[i+1] = y[i] + f(tempos[i])*delta_t
        tempos[i+1] = tempos[i] + delta_t

    erro = abs(solucao(1) - y[n-1])
    
    return delta_t, erro, tempos, y

#-----------------------------------------------------------------------------

def solucao(t):
    '''Solucaoo da EDO calculada na questao 1 da tarefa 1'''
    s = math.exp(t)*math.sin(2*t) + 1
    return s


#############################################################################
#############################################################################
    

def main():
    ''' Programa principal que realiza o metodo de Euler k vezes onde a cada
    realizacao o intervalo utilizado no metodo e da forma 1/2^k'''
    k = 14 #Quantidade de iteracoes
    deltas = [0]*(k-2)
    erros = [0]*(k-2)
    for i in range(2,k):
         n = 2**i
         metodo = euler(0,1,n)
         deltas[i-2], erros[i-2] = metodo[0:2]
         plt.plot(metodo[2],metodo[3], color='black')
         plt.xlabel('t')
         plt.ylabel('y(t)')
         plt.title('Convergência do Método de Euler')
         



    print(" Passo(Delta t) | Erro Absoluto")
    print("="*39)
    
    for i in range(len(deltas)):
        print("    %f    |"%deltas[i], end="")
        print("    %f  "%erros[i])
    print("="*39)
    
    plt.show()
    
    X = np.linspace(0, 1, 1024)
    Y = np.linspace(0, 1, 1024)
    for i in range(len(X)):
        Y[i] = solucao(X[i])
    aprox = euler(0,1,n)
    
    plt.plot(aprox[2],aprox[3], color='black', linestyle=(0,(1,1,3,1)),
             label = 'Aproximação numérica com n = 8196')
    plt.plot(X, Y, c = 'k', label = "Função verdadeira conhecida")
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title('Comparação entre o Método Numérico e a Solução Exata')
    plt.legend()
    plt.show()

main()
