#Phase Space representation of th bicycle model#
#The idea of this is to giev na intative look into te behaviour of the bicycle model under certian conditions by using a 'phase space like' representation
# not entirley sre how this will work atm but will bash into it and see what comes of it.
# https://thef1clan.com/2020/12/23/vehicle-dynamics-the-dynamic-bicycle-model/

import numpy as np
import matplotlib.pyplot as plt

Iz = 75; #yaw inertia
a = 0.81355; #front wheel moment arm
b = 0.7125; # rear wheel moment arm
Cs = 300; # Cornering stiffness n/deg
m = 300;  #mass kg
tw = 1.2 #track width


Sm = np.zeros((2,2)); #system
u = np.zeros((2,1)); #input matrix


vel = 10; #velocity
drag = 0.5*1.001*1.01*1.6*vel**2 #drag
Delta = 0; #Steer angle


#initilaising the system matrix of the state space representaion
#Sm[1,1] = (2*Cs)/(m*vel);
#Sm[1,2] = (Cs*a - Cs*b)/(m*vel) + vel;
#Sm[1,2] = (Cs*a - Cs*b)/(vel*Iz);
#Sm[2,2] = (Cs*a**2 + Cs*b**2)/(vel*Iz);

#Sm = -1*Sm;

#Initialising the input matrix
#u[1,1] = Cs/m;
#u[1,2] = (Cs*a)/Iz;
step = 0.1; 
betamax = 20; # body slip
yawratemax = 100; # yaw rate

beta = np.linspace(-betamax,betamax,10)
yawrate = np.linspace(-yawratemax,yawratemax,10)
alpha = np.zeros((len(beta),len(yawrate)))
BetaDot = np.zeros((len(beta),len(yawrate)))

print(beta)

for i in range(len(beta)): 
  for j in range(len(yawrate)):

   Beta = beta[i]
   Yawr = yawrate[j]

   alphaf = Beta + ((Yawr*a)/vel) - Delta # front slip angle
   alphar = Beta - ((Yawr*b)/vel) # rear slip angle

   fyf = 2*alphaf*Cs
   fyr = 2*alphar*Cs
   fxr = drag/2

   BetaDot[i,j] = ((fyf+fyr)/(m*vel)) - Yawr
   alpha[i,j] = (fyf*a - fyr*b)/Iz

   
print(BetaDot)

plt.quiver(beta, yawrate, BetaDot, alpha)
plt.xlabel('body slip [deg]')
plt.ylabel('Yaw Rate [deg/s]')
plt.title('Bicycle Model Phase Space')
plt.show()



  



