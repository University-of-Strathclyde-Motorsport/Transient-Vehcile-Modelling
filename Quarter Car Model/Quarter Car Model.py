import numpy as np
import matplotlib.pyplot as plt

dt = 0.001;

tfin = 5;

time = np.arange(0,tfin+dt,dt); #initialise time array
Ss = np.zeros(len(time)+1); #sprung displacmenet
Vs = np.zeros(len(time)+1);
As = np.zeros(len(time));
Sus = np.zeros(len(time)+1);
Vus = np.zeros(len(time)+1);
Aus = np.zeros(len(time));

Vs[0] = -4.43; #initialis first velocity value
Vus[0] = -4.43; #initialis first velocity value

#initilaise spring and damper values

Ms = 61.1; #sprung mass
Mus = 10; # unsprung mass
Keff = 30000; # effective spring rate
Cs = 1000; # dmaping coefficent
Kt = 100000; # tyre rate 


# calculation loop
for i in range(len(time)):

  Aus[i] = (-Kt*Sus[i] - Keff*(Sus[i]-Ss[i]) - Cs*(Vus[i] - Vs[i]))/Mus;
  As[i] = (-Keff*(Ss[i]-Sus[i]) - Cs*(Vs[i]-Vus[i]))/Ms;

  Vus[i+1] = Vus[i] + Aus[i]*dt;
  Vs[i+1] = Vs[i] + As[i]*dt;

  Sus[i+1] = Sus[i] + Vus[i]*dt;
  Ss[i+1] = Ss[i] + Vs[i]*dt;


plt.plot(time,Sus[0:len(time)],label = "Unsprung Mass Displacement")
plt.plot(time,Ss[0:len(time)],label = "Sprung Mass Displacement")
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title('Diaplacement Plot')
plt.show()





