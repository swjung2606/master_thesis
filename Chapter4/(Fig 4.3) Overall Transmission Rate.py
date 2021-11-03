from math import *
import numpy as np
import matplotlib.pyplot as plt


t=np.linspace(0,478,479)
h=t[1]

aa=[]
R0=[]

a1=1.2701
a2=0.1948
a3=0.00011
k1=0.7139
k2=0.0752
b1=0.1417
b2=-2.2807


c1=0.3847
c2=0.2961
l1=60

for i in range(0,187):
    aa.append( ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t[i])) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t[i]+b2)) ) )

for i in range(187,479):
    aa.append( ( 0.5*(c1-c2)*cos( (t[i]-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) ) )



plt.plot(t,aa,color="blue")
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.ylabel("Transmission Rate $\u03B2(t)$")
plt.title("Transmission Rate")
plt.legend(["Overall Transmission Rate"])
plt.annotate("Social Distancing \nCampaign By \nCentral Government \nat t=187",xy=(187,0.3228),xytext=(187,0.8),arrowprops=dict(facecolor="black",shrink=0.05))
plt.annotate("Ease the Campaign By Central Government at t=237",xy=(237,0.3185),xytext=(237,0.12),arrowprops=dict(facecolor="black",shrink=0.05))
plt.annotate("Social Distancing \nCampaign By \nCentral Government \nat t=287",xy=(287,0.3541),xytext=(287,1.12),arrowprops=dict(facecolor="black",shrink=0.05))
plt.grid(True)
plt.show()
