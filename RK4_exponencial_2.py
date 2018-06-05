import numpy as np #importa el paquete que permite usar funciones matemรกticas
import pylab as plb #importa pquete que nos permite hacer grรกficos

def int_rk4(f,x,dt,params):
    
    k_1 = f(x,params)
    k_2 = f(x+dt*0.5*k_1,params)
    k_3 = f(x+dt*0.5*k_2,params)
    k_4 = f(x+dt*k_3,params)
    y=x + dt*(k_1/6 + k_2/3 + k_3/3 + k_4/6)
    return y

#quiero integrar la funciรณn dx/dt=ax, para eso me construyo la funciรณn derivada que le voy a dar al RK4
def derivada(x,params):

    a=params #asรญ lo elijo de afuera
    dx=a*x   #defino la derivada de la exponencial

    return np.array(dx)  #le pido a la funciรณn que me devuelva la derivada  

#ahora integro para el tiempo que yo quiera

#me armo el vector de tiempo
t0=0 #tiempo inicial
tf=8.1 #tiempo final
dt=0.1 #paso del tiempo
t=np.arange(t0,tf,dt) #vector que va de t0 a tf con un paso dt

#le doy un valor a mi parรกmetro
a=np.log(2) #El 2 significa que por cada unidad de tiempo la poblaciรณn se duplica         

#me armo la variable que voy a integrar
           
x=[] # x es un vector de ceros que vamos a ir llenando
x.append(1) #fijo condiciรณn inicial

#integro usando el mรฉtodo RK4 definido antes

for i in range (len(t)-1): #"para todos los i es menores a la longitud de t"
    
    x.append(int_rk4(derivada,x[i],dt,a)) #le doy el tiempo y la x correspondientes al paso i-esimo para calcular el paso siguiente
    
#luego del for tenemos la ecuaciรณn integrada para t entre 0 y tf, comparemos con la soluciรณn teรณrica

xteo=np.exp(a*t)

#ploteemos las funciones a ver si dan lo mismo

plb.plot(t,xteo,'-',lw=5,label='Teórica')  #la rayita '-' hace que grafique con una raya, lw es el grosor de la línea y el label la leyenda
plb.plot(t,x,'.',label='RK4')   #el puntito '.' hace que grafique con puntitos             
plb.xlabel('tiempo')
plb.ylabel('Población')
plb.legend(loc='upper left',fontsize=10) #le digo dรณnde poner las leyendas
plb.show()

