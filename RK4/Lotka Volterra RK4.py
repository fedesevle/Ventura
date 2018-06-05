import numpy as np #importa el paquete que permite usar funciones matemáticas
import pylab as plb #importa pquete que nos permite hacer gráficos

def int_rk4(f,x,dt,params):              #la misma función de RK4 me sirve para integrar muchas variables, ahora x va a ser un vector
    
    k_1 = f(x,params)
    k_2 = f(x+dt*0.5*k_1,params)
    k_3 = f(x+dt*0.5*k_2,params)
    k_4 = f(x+dt*k_3,params)
    y=x + dt*(k_1/6 + k_2/3 + k_3/3 + k_4/6)
    return y

#quiero integrar la función dx/dt=ax, para eso me construyo la función derivada que le voy a dar al RK4
def derivada(x,params):

    a=params[0]
    b=params[1]
    c=params[2]
    d=params[3]
    
    x0=x[0] #primer especie
    x1=x[1] #segunda especie
    
    dx0=a*x0-b*x0*x1   #ecuación diferencial para la primer especie, presa
    dx1=-c*x1+d*x0*x1   #ecuación diferencial para la segunda especie, predador

    return np.array([dx0,dx1])  #le pido a la función que me devuelva la derivada
    
    
    #defino parámetros

a=1
b=0.01
c=1
d=0.01

params=[a,b,c,d] #armos vector de parámetros

#armos vector de tiempo

t0=0
tf=50
dt=0.001
t=np.arange(t0,tf,dt)

#armo la x que va a contener la información de las dos especies

x=np.zeros([2,len(t)]) #ahora quiero que sea una matriz de 2xlen(t)

delta=0.1 # Oscilaciones sinusoidales
#delta=0.75 # Oscilaciones deformadas
#delta=0.9999 # Oscilaciones criticas
#delta=1 # Extinsion de los depredadores
p#delta=-1 # Extinsion de las presas y luego de los depredadores
x[:,0]=[a/b*(1+delta),c/d*(1-delta)] #fijo condiciones iniciales


  
#ya podemos integrar de forma muy parecida a una variable

for i in range (len(t)-1):
        x[:,i+1]=int_rk4(derivada,x[:,i],dt,params) #ahora donde ponía números pongo vectorcitos en las x
    
#veo como dio

x0=x[0,:] #resultados de la primer especio
x1=x[1,:] #resultados de la segunda especie

#grafiquemos 
plb.subplot(211)   
plb.plot(t,x0,label='Presas')  #le pongo leyendas a cada especie    
plb.plot(t,x1,label='Deredadores')     
plb.xlabel('Tiempo')
plb.ylabel('Población')
plb.title('Modelo Lotka-Volterra')
plb.legend(loc='upper right',fontsize=10) #le digo dónde poner las leyendas
plb.show()  #le digo que muestre las leyendas (sin este comando no van a aparecer)      plb.subplot(121)   
plb.subplot(212)   
plb.plot(x0,x1)  #le pongo leyendas a cada especie  
plb.xlabel('Presas')
plb.ylabel('Depredadores')
