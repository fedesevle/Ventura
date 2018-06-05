import numpy as np #importa el paquete que permite usar funciones matemáticas
import pylab as plb #importa pquete que nos permite hacer gráficos

a=1
b=0.01
c=1
d=0.01
params=[a,b,c,d] #armo vector de parámetros

x=[]
y=[]
delta=0 # Oscilaciones sinusoidales
#delta=0.1 # Oscilaciones sinusoidales
#delta=0.75 # Oscilaciones deformadas
#delta=1 # Extinsion de los depredadores
#delta=-1 # Extinsion de las presas y luego de los depredadores
x.append(a/b*(1+delta))
y.append(c/d*(1-delta))
t=[]
t.append(0)
tf=50

def Gillespie_LV(x,y,t,tf,params):
    i=0 
    a=params[0]
    b=params[1] 
    c=params[2] 
    d=params[3]       
    while t[len(t)-1]<tf:
        dR_xx=a*x[i]
        dR_yy=c*y[i]
        dR_xy=b*x[i]*y[i]
        dR_yx=d*x[i]*y[i]
        dR=np.sum(dR_xx+dR_xy+dR_yx+dR_yy)
        dt=np.random.exponential(1/dR)
        p=np.random.rand()
        if 0<p and p<dR_xx/dR:
            x.append(x[i]+1)
            y.append(y[i])
        elif dR_xx/dR<p and p<(dR_xx+dR_xy)/dR:
            x.append(x[i]-1)
            y.append(y[i])
        elif (dR_xx+dR_xy)/dR<p and p<(dR_xx+dR_xy+dR_yx)/dR:
            x.append(x[i])
            y.append(y[i]+1)
        else:
            x.append(x[i])
            y.append(y[i]-1)              
        t.append(t[i]+dt)
        i=i+1
    return x,y,t

x,y,t=Gillespie_LV(x,y,t,tf,params)

plb.subplot(211)
plb.plot(t,x,label='Presa')
plb.plot(t,y,label='Depredador')
plb.xlabel('Tiempo')
plb.ylabel('Población')
plb.title('Modelo Lotka-Volterra')
plb.legend(loc='upper left',fontsize=10)  #le digo dónde poner las leyendas
plb.show()  #le digo que muestre las leyendas (sin este comando no van a aparecer)
plb.subplot(212)
plb.plot(x,y)
plb.xlabel('Presas')
plb.ylabel('Depredadores')
