import numpy as np
import matplotlib.pyplot as plt

def eulerExplicito(n,h,i0) :
  i = np.empty((n+1))
  i[0] = i0
  tempos = np.empty((n+1))
  tempos[0] = 0

  for j in range(1,n+1):
    i[j] = (1 + h*(-5)) * i[j-1] #Euler explícito
    tempos[j] = tempos[j-1] + h



  return tempos, i


def main():
  k= 5 #quantidade de iterações
  t0 = 0 #tempo inicial
  tf = 1 #tempo final
  i0 = 2 #corrente inicial
  erroAnterior = 0 
  n =8 #numero de passos inicial
  cores = ['dimgray','gray','darkgray','darkgrey','silver']
  resultadoexato = 2*np.exp(-5*tf) # di(t)/dt = -5i(t) ; i(0)=2  --->   i(t) = 2e^(-5t)
  for j in range(1,k+1):
    h = (tf-t0)/n  
    resultadoEuler = eulerExplicito(n,h,i0)
    plt.plot(resultadoEuler[0],resultadoEuler[1], color=cores[j-1])
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title('Convergência do Método de Euler')
    erro = np.absolute(resultadoexato - resultadoEuler[1][-1])

    if(j == 1):
      print(n," & ", h," & ",erro," & ","---------------" )
    

    else:
      ordem = np.absolute(np.log2(np.absolute(erro / erroAnterior)))
      print(n," & ", h," & ",erro," & ",ordem )   
    
    n *= 2
    erroAnterior = erro
    
  plt.show()    
main()