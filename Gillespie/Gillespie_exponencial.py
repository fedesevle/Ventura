import numpy as np #importa el paquete que permite usar funciones matemรกticas
import pylab as plb #importa pquete que nos permite hacer grรกficos

# Defino el parametro
a=np.log(2) #El 2 significa que por cada unidad de tiempo la poblaciรณn se duplica   

x=[] # x es el vector poblacion que iremos llenando
x.append(1) # Poblacion inicial

t=[] # t es el tiempo que iremos llenando
t.append(0) # Tiempo inicial
tf=8 # tiempo final hasta donde correrá Gillespie

# Defino la función que es el algoritmo de Gillespie
def Gillespie(x,t,tf):     
    i=0        
    while t[len(t)-1]<tf:
        dt=np.random.exponential(1/(a*x[i]))
        dt=1/(a*(x[i]+1/2))
        x.append(x[i]+1)
        t.append(t[i]+dt)
        i=i+1
    return x,t
    
# Aplico la función
x,t=Gillespie(x,t,tf)
    
# Calculo la solución teórica - determinista
xteo=x[0]*np.exp(a*np.asarray(t))

# Graficamos la población teórica y la obtenida
#plb.plot(t,xteo,'-',linewidth=5,label='Determinista')  #la rayita al final hace que grafique con una raya 
plb.plot(t,x,'.',linewidth=5,label='Estocástico')   #el puntito al final hace que grafique con puntitos           
plb.xlabel('tiempo')
plb.ylabel('Población')
plb.legend(loc='upper left',fontsize=10)
plb.show()
