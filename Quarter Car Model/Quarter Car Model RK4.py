## RK4 version ##
#quarter car model solved using RK4 instead of gauss forward method
import numpy as np
import matplotlib.pyplot as plt

dt = 0.001;

tfin = 5;

time = np.arange(0,tfin,dt); #initialise time array
Ss = np.zeros(len(time)); #sprung displacmenet
Vs = np.zeros(len(time));
As = np.zeros(len(time));
Sus = np.zeros(len(time));
Vus = np.zeros(len(time));
Aus = np.zeros(len(time));

Vs[0] = -4.43; #initialis first velocity value
Vus[0] = -4.43; #initialis first velocity value

#initilaise spring and damper values

Ms = 61.1; #sprung mass
Mus = 10; # unsprung mass
Keff = 30000; # effective spring rate
Cs = 1000; # dmaping coefficent
Kt = 100000; # tyre rate 

#define functions that are used in the RK4 loop

def Function_f1(Vus):
    #dx1/dt = v1
   return Vus  
   

def Function_f2(Kt,Keff,Cs,Mus,Sus,Ss,Vus,Vs):
    #dv1/dt equatiom
    return (-Kt*Sus - Keff*(Sus-Ss) - Cs*(Vus - Vs))/(Mus)
 

def Function_f3(Vs):
    return Vs
  

def Function_f4(Keff,Cs,Ms,Sus,Ss,Vus,Vs):
    return (-Keff*(Ss-Sus) - Cs*(Vs - Vus))/Ms
    

#now that we have defined the RK4 functions we can use them in the RK4 Calculation loop

for i in range(len(time)-1):
    k1 = dt*Function_f1(Vus[i])
    l1 = dt*Function_f2(Kt,Keff,Cs,Mus,Sus[i],Ss[i],Vus[i],Vs[i])
    m1 = dt*Function_f3(Vs[i])
    n1 = dt*Function_f4(Keff,Cs,Ms,Sus[i],Ss[i],Vus[i],Vs[i])
    
    k2 = dt*Function_f1(Vus[i] + k1*0.5)
    l2 = dt*Function_f2(Kt , Keff , Cs , Mus , Sus[i] + k1*0.5 , Ss[i] + 0.5*m1 , Vus[i] + 0.5*l1 , Vs[i] + 0.5*n1)
    m2 = dt*Function_f3(Vs[i])
    n2 = dt*Function_f4(Keff,Cs,Ms,Sus[i] + 0.5*k1 , Ss[i] + 0.5*m1 , Vus[i] + 0.5*l1 , Vs[i] + 0.5*n1)
    
    k3 = dt*Function_f1(Vus[i])
    l3 = dt*Function_f2(Kt,Keff,Cs,Mus,Sus[i] + 0.5*k2 , Ss[i] + 0.5*m2 , Vus[i] + 0.5*l2 , Vs[i] + 0.5*n2)
    m3 = dt*Function_f3(Vs[i])
    n3 = dt*Function_f4(Keff,Cs,Ms,Sus[i] + 0.5*k2 , Ss[i] + 0.5*m2 , Vus[i] + 0.5*l2 , Vs[i] + 0.5*n2)
    
    k4 = dt*Function_f1(Vus[i])
    l4 = dt*Function_f2(Kt , Keff , Cs , Mus , Sus[i] + k3 , Ss[i] + m1 , Vus[i] + l3 , Vs[i] + n3)
    m4= dt*Function_f3(Vs[i])
    n4 = dt*Function_f4(Keff , Cs , Ms , Sus[i] + k3 , Ss[i] + m1 , Vus[i] + l3 , Vs[i] + n3)
    
    Sus[i+1] = Sus[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    Vus[i+1] = Vus[i] + (1/6)*(l1 + 2*l2 + 2*l3 + l4)
    Ss[i+1] = Ss[i] + (1/6)*(m1 + 2*m2 + 2*m3 + m4)
    Vs[i+1] = Vs[i] + (1/6)*(n1 + 2*n2 + 2*n3 + n4)
    
    del k1, k2, k3 ,k4 ,l1 , l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4
    



plt.plot(time,Sus)
plt.plot(time,Ss)
plt.title('RK4 Method')
plt.show()
