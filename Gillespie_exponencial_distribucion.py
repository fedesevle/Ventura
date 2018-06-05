import numpy as np #importa el paquete que permite usar funciones matemรกticas
import pylab as plb #importa pquete que nos permite hacer grรกficos

#Le doy un valor a mi parametro
a=np.log(2) #El 2 significa que por cada unidad de tiempo la poblacion se duplica   

# Defino la función de Gillespie hasta que alcanza el tiempo final tf
def Gillespie(x,t,tf):     
    i=0        
    while t[len(t)-1]<tf:
        dt=np.random.exponential(1/(a*x[i]))
        x.append(x[i]+1)
        t.append(t[i]+dt)
        i=i+1
    return x,t

# Calculo muchos crecimientos de la población para obtener la distribución
X=[] #acá voy a guardar las poblaciones
T=[] #acá voy a guardar los tiempos
N=2000 #Cantidad de repeticiones

for i in range(N):
    x=[] # x es el número de individuos que van a ir creciendo con el tiempo
    x.append(1) # x inicial (1 individuo)
    t=[] # t es el tiempo
    t.append(0) # tiempo inicial (0)
    tf=8 # tiempo final hasta donde crece la población
    x,t=Gillespie(x,t,tf)
    X.append(x)
    T.append(t)


#Veamos si la distribución de las poblaciones en función del tiempo se parece a la curva teórica

# Interpolo todas las poblaciones al mismo tiempo, para compararlas
X_interp=[]
T_interp=np.arange(0,tf,0.1)
for i in range(N):
    X_interp.append(np.interp(T_interp,T[i],X[i]))
X_matrix=np.array(X_interp)

# Calculo las medias y los quantiles 5, 25, 75 y 95 para graficarlos
X_mean=np.zeros(len(T_interp))
X_up=np.zeros(len(T_interp))
X_down=np.zeros(len(T_interp))
X_upp=np.zeros(len(T_interp))
X_downn=np.zeros(len(T_interp))
for i in range(len(T_interp)):
    X_mean[i]=np.mean(X_matrix[:,i])
    X_up[i]=np.percentile(X_matrix[:,i],75)
    X_down[i]=np.percentile(X_matrix[:,i],25)
    X_upp[i]=np.percentile(X_matrix[:,i],97.5)
    X_downn[i]=np.percentile(X_matrix[:,i],2.5)
    
#Curva teórica
xteo=x[0]*np.exp(a*T_interp)

# Grafico
plb.plot(T_interp,X_mean,'-',label='Media',linewidth=5)   #el puntito al final hace que grafique con puntitos     
plb.show()      
plb.fill_between(T_interp,X_down,X_up,alpha=0.25,label='Intervalo de Confianza 50%')   #el puntito al final hace que grafique con puntitos
plb.fill_between(T_interp,X_downn,X_upp,alpha=0.10,label='Intervalo de Confianza 95%')   #el puntito al final hace que grafique con puntitos 
plb.plot(T_interp,xteo,'-.',linewidth=3,label='Determinista')  #la rayita al final hace que grafique con una raya  
plb.xlabel('tiempo')
plb.ylabel('N')
for i in range(50):
    plb.plot(T[i],X[i],alpha=0.25)
plb.plot(T[i],X[i],alpha=0.25,label='Estocásticos')
plb.xlim([0,tf])        
plb.ylim([0,350])     
plb.legend(loc='upper left',fontsize=10)
plb.show()

