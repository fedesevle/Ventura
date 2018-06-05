# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:44:04 2018

@author: USUARIO
"""


import numpy as np #importa el paquete que permite usar funciones matemáticas
import pylab as plb #importa pquete que nos permite hacer gráficos

def int_rk4(f,x,dt,params):
    
    k_1 = f(x,params)
    k_2 = f(x+dt*0.5*k_1,params)
    k_3 = f(x+dt*0.5*k_2,params)
    k_4 = f(x+dt*k_3,params)
    y=x + dt*(k_1/6 + k_2/3 + k_3/3 + k_4/6)
    return y

#quiero integrar la función dx/dt=ax, para eso me construyo la función derivada que le voy a dar al RK4
def derivada(x,params):

    a=params #así lo elijo de afuera
    dx=a*x   #defino la derivada de la exponencial

    return np.array(dx)  #le pido a la función que me devuelva la derivada  

#ahora integro para el tiempo que yo quiera

#me armo el vector de tiempo
t0=0 #tiempo inicial
tf=8 #tiempo final
dt=0.1 #paso del tiempo
t=np.arange(t0,tf,dt) #vector que va de t0 a tf con un paso dt

#le doy un valor a mi parámetro
a=np.log(2) #Significa que por cada unidad de tiempo la población se duplica         

#me armo la variable que voy a integrar
           
x=np.zeros(len(t)) # x es un vector de ceros que vamos a ir llenando
x[0]=1 #fijo condición inicial

#integro usando el método RK4 definido antes

for i in range (len(t)-1): #"para todos los i es menores a la longitud de t"
    
    x[i+1]=int_rk4(derivada,x[i],dt,a) #le doy el tiempo y la x correspondientes al paso i-esimo para calcular el paso siguiente
    
#luego del for tenemos la ecuación integrada para t entre 0 y tf, comparemos con la solución teórica

xteo=np.exp(a*t)

#ploteemos las funciones a ver si dan lo mismo

plb.plot(t,x,'.')   #el puntito al final hace que grafique con puntitos           
plb.plot(t,xteo,'-')  #la rayita al final hace que grafique con una raya  
plb.xlabel('tiempo')
plb.ylabel('N')
        
          
               
