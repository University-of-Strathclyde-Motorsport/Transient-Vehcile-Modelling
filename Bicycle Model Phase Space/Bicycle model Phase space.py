#Phase Space representation of th bicycle model#
#The idea of this is to giev na intative look into te behaviour of the bicycle model under certian conditions by using a 'phase space like' representation
# not entirley sre how this will work atm but will bash into it and see what comes of it.
# RCVD, chapter 5

import numpy as np
import matplotlib.pyplot as plt

Iz = 75; #yaw inertia
a = 0.81355; #front wheel moment arm
b = 0.7125; # rear wheel moment arm
Csf = 300; # Cornering stiffness n/deg
Csr = 350; # rear cornering stiffness n/deg
m = 300;  #mass kg


#tw = 1.2 #track width
#Sm = np.zeros((2,2)); #system
#u = np.zeros((2,1)); #input matrix


vel = 15; #velocity
drag = 0.5*1.001*1.01*1.6*vel**2 #drag
Delta = 10; #Steer angle

betamax = 20; # body slip deg
yawratemax = 100; # yaw rate deg/s

beta = np.linspace(-betamax,betamax,10)
yawrate = np.linspace(-yawratemax,yawratemax,10)
alpha = np.zeros((len(beta),len(yawrate)))
BetaDot = np.zeros((len(beta),len(yawrate)))


for i in range(len(beta)): 
  for j in range(len(yawrate)):

   Beta = beta[i]
   Yawr = yawrate[j]

   alphaf = Beta + ((Yawr*a)/vel) - Delta # front slip angle
   alphar = Beta - ((Yawr*b)/vel) # rear slip angle

   fyf = 2*alphaf*Csf
   fyr = 2*alphar*Csr
   fxr = drag/2
   
   ## add in bilinear tyre model
   #the tyre willa ct as linear up to 8 deg, then it will level off
   if abs(fyf) > Csf*8:
     if fyf > 0:
       fyf = Csf*8
     else:
       fyf = -Csf*8
   
   
   if abs(fyr) > Csr*8:
     if fyf > 0:
       fyf = Csr*8
     else:
       fyf = -Csr*8
     
     
    
   BetaDot[i,j] = ((fyf+fyr)/(m*vel)) - Yawr
   alpha[i,j] = (fyf*a - fyr*b)/Iz





plt.quiver(beta, yawrate, BetaDot, alpha)
plt.xlabel('body slip [deg]')
plt.ylabel('Yaw Rate [deg/s]')
plt.title('Bicycle Model Phase Space')
plt.show()

## conculsion ##
#algo works but needs MATLAB for Tyre model, i'm not building one in python
#still cool none the less,
#This is basically the bases of a YMD script



  



