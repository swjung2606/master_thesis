# Logistics model sensitivities with respect to beta


from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


Sy=[31451865]
Sm=[14578212]
So=[5194184]
Ey=[217]
Em=[577]
Eo=[30]
A=[0]
Iy=[0]
Im=[0]
Io=[0]
My=[18]
Mm=[11]
Mo=[2]
Dy=[0]
Dm=[0]
Do=[0]
R=[0]
IT=[31]


S_Sy=[0]
S_Sm=[0]
S_So=[0]
S_Ey=[0]
S_Em=[0]
S_Eo=[0]
S_A=[0]
S_Iy=[0]
S_Im=[0]
S_Io=[0]
S_My=[0]
S_Mm=[0]
S_Mo=[0]
S_Dy=[0]
S_Dm=[0]
S_Do=[0]
S_R=[0]
S_IT=[0]



S_Sy_x=[0]
S_Sm_x=[0]
S_So_x=[0]
S_Ey_x=[0]
S_Em_x=[0]
S_Eo_x=[0]
S_A_x=[0]
S_Iy_x=[0]
S_Im_x=[0]
S_Io_x=[0]
S_My_x=[0]
S_Mm_x=[0]
S_Mo_x=[0]
S_Dy_x=[0]
S_Dm_x=[0]
S_Do_x=[0]
S_R_x=[0]
S_IT_x=[0]



N = 51225116

q = 0.01
eta = 1/5.1
p = 0.844
gamma = 1/14
alpha = 1/17.5
x = 1/4 
x1 = 1/3
dy = 0.0006 # Based on 8/22 number of death
dm = 0.0104
do = 0.1341
theta = 1/13
epsilon = 1.1254



a1=1.2701
a2=0.1948
a3=0.00011
k1=0.7139
k2=0.0752
b1=0.1417
b2=-2.2807


# Solution equations part

def Sy_prime(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Sm_prime(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def So_prime(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Ey_prime(t,Sy,Ey,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey

def Em_prime(t,Sm,Em,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em

def Eo_prime(t,So,Eo,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo

def A_prime(t,Ey,Em,Eo,A):
        return (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A

def IT_prime(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def Iy_prime(t,Ey,Iy):
        return (eta)*p*Ey - x*Iy

def Im_prime(t,Em,Im):
        return (eta)*p*Em - x*Im

def Io_prime(t,Eo,Io):
        return (eta)*p*Eo - x*Io

def My_prime(t,Iy,My):
        return  x*Iy - ( alpha*(1-dy)+dy*theta )*My

def Mm_prime(t,Im,Mm):
        return  x*Im - ( alpha*(1-dm)+dm*theta )*Mm

def Mo_prime(t,Io,Mo):
        return  x*Io - ( alpha*(1-do)+do*theta )*Mo

def Dy_prime(t,My):
        return dy*(theta)*My

def Dm_prime(t,Mm):
        return dm*(theta)*Mm

def Do_prime(t,Mo):
        return do*(theta)*Mo

def R_prime(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )


# Sensitivity equations part for beta

def Sy_1(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sy,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sy/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sy*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) -(Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def Sm_1(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sm,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sm/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sm*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) -(Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def So_1(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_So,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -(epsilon)*( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_So/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (So*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) -(epsilon)*(So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def Ey_1(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sy,S_Ey,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sy/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sy*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Ey + (Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def Em_1(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sm,S_Em,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sm/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sm*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Em + (Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def Eo_1(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_So,S_Eo,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return (epsilon)*( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_So/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (So*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Eo + (epsilon)*(So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(Iy+Im+Io+q*A)

def A_1(t,S_Ey,S_Em,S_Eo,S_A):
        return (eta)*(1-p)*(S_Ey+S_Em+S_Eo)-(gamma)*S_A

def Iy_1(t,S_Ey,S_Iy):
        return (eta)*p*S_Ey - x*S_Iy

def Im_1(t,S_Em,S_Im):
        return (eta)*p*S_Em - x*S_Im

def Io_1(t,S_Eo,S_Io):
        return (eta)*p*S_Eo - x*S_Io

def My_1(t,S_Iy,S_My):
        return x*S_Iy - ((alpha)*(1-dy) + dy*(theta))*S_My

def Mm_1(t,S_Im,S_Mm):
        return x*S_Im - ((alpha)*(1-dm) + dm*(theta))*S_Mm

def Mo_1(t,S_Io,S_Mo):
        return x*S_Io - ((alpha)*(1-do) + do*(theta))*S_Mo

def Dy_1(t,S_My):
        return dy*(theta)*S_My

def Dm_1(t,S_Mm):
        return dm*(theta)*S_Mm

def Do_1(t,S_Mo):
        return do*(theta)*S_Mo

def R_1(t,S_A,S_My,S_Mm,S_Mo):
        return (gamma)*S_A + (alpha)*( (1-dy)*S_My + (1-dm)*S_Mm + (1-do)*S_Mo )

def IT_1(t,S_Ey,S_Em,S_Eo):
        return (eta)*p*(S_Ey+S_Em+S_Eo)



# Sensitivity equations part for x

def Sy_2(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sy,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sy/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sy*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) 

def Sm_2(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sm,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sm/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sm*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) 

def So_2(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_So,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return -(epsilon)*( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_So/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (So*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) 

def Ey_2(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sy,S_Ey,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sy/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sy/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sy*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Ey 

def Em_2(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_Sm,S_Em,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_Sm/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (Sm/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (Sm*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Em 

def Eo_2(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do,S_So,S_Eo,S_A,S_Iy,S_Im,S_Io,S_My,S_Mm,S_Mo,S_Dy,S_Dm,S_Do):
        return (epsilon)*( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*( (Iy+Im+Io+q*A)*(S_So/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So*q*S_A/(N-(My+Mm+Mo+Dy+Dm+Do))) + (So/(N-(My+Mm+Mo+Dy+Dm+Do)))*(S_Iy+S_Im+S_Io) + (So*(Iy+Im+Io+q*A))*(S_My+S_Mm+S_Mo+S_Dy+S_Dm+S_Do)/((N-(My+Mm+Mo+Dy+Dm+Do))**2) ) - (eta)*S_Eo

def A_2(t,S_Ey,S_Em,S_Eo,S_A):
        return (eta)*(1-p)*(S_Ey+S_Em+S_Eo)-(gamma)*S_A

def Iy_2(t,Iy,S_Ey,S_Iy):
        return (eta)*p*S_Ey - x*S_Iy -Iy

def Im_2(t,Im,S_Em,S_Im):
        return (eta)*p*S_Em - x*S_Im -Im

def Io_2(t,Io,S_Eo,S_Io):
        return (eta)*p*S_Eo - x*S_Io -Io

def My_2(t,Iy,S_Iy,S_My):
        return x*S_Iy - ((alpha)*(1-dy) + dy*(theta))*S_My +Iy

def Mm_2(t,Im,S_Im,S_Mm):
        return x*S_Im - ((alpha)*(1-dm) + dm*(theta))*S_Mm +Im

def Mo_2(t,Io,S_Io,S_Mo):
        return x*S_Io - ((alpha)*(1-do) + do*(theta))*S_Mo +Io

def Dy_2(t,S_My):
        return dy*(theta)*S_My

def Dm_2(t,S_Mm):
        return dm*(theta)*S_Mm

def Do_2(t,S_Mo):
        return do*(theta)*S_Mo

def R_2(t,S_A,S_My,S_Mm,S_Mo):
        return (gamma)*S_A + (alpha)*( (1-dy)*S_My + (1-dm)*S_Mm + (1-do)*S_Mo )

def IT_2(t,S_Ey,S_Em,S_Eo):
        return (eta)*p*(S_Ey+S_Em+S_Eo)




t=np.linspace(0,186,187)
h=t[1]


for i in range(0,186):  
    Sy.append( Sy[i] + h*Sy_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Sm.append( Sm[i] + h*Sm_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    So.append( So[i] + h*So_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Ey.append( Ey[i] + h*Ey_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Em.append( Em[i] + h*Em_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Eo.append( Eo[i] + h*Eo_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    A.append( A[i] + h*A_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]) ) )
    Iy.append( Iy[i] + h*Iy_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]) ) )
    Im.append( Im[i] + h*Im_prime(t[i]+h/2, Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]) ) )
    Io.append( Io[i] + h*Io_prime(t[i]+h/2, Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]) ) )
    My.append( My[i] + h*My_prime(t[i]+h/2, Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]) ) )
    Mm.append( Mm[i] + h*Mm_prime(t[i]+h/2, Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]) ) )
    Mo.append( Mo[i] + h*Mo_prime(t[i]+h/2, Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]) ) )
    Dy.append( Dy[i] + h*Dy_prime(t[i]+h/2, My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]) ) )
    Dm.append( Dm[i] + h*Dm_prime(t[i]+h/2, Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]) ) )
    Do.append( Do[i] + h*Do_prime(t[i]+h/2, Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]) ) )   
    R.append( R[i] + h*R_prime(t[i]+h/2, A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]) ) )
    IT.append( IT[i] + h*IT_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )


    S_Sy.append( S_Sy[i] + h*Sy_1(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sy[i]+(h/2)*Sy_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_Sm.append( S_Sm[i] + h*Sm_1(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sm[i]+(h/2)*Sm_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_So.append( S_So[i] + h*So_1(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_So[i]+(h/2)*So_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_Ey.append( S_Ey[i] + h*Ey_1(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sy[i]+(h/2)*Sy_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Ey[i]+(h/2)*Ey_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_Ey[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_Em.append( S_Em[i] + h*Em_1(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sm[i]+(h/2)*Sm_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Em[i]+(h/2)*Em_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_Em[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_Eo.append( S_Eo[i] + h*Eo_1(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_So[i]+(h/2)*So_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Eo[i]+(h/2)*Eo_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_Eo[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]), S_Dy[i]+(h/2)*Dy_1(t[i],S_My[i]), S_Dm[i]+(h/2)*Dm_1(t[i],S_Mm[i]), S_Do[i]+(h/2)*Do_1(t[i],S_Mo[i]) ) )
    S_A.append( S_A[i] + h*A_1(t[i]+h/2, S_Ey[i]+(h/2)*Ey_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_Ey[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Em[i]+(h/2)*Em_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_Em[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Eo[i]+(h/2)*Eo_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_Eo[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]) ) )
    S_Iy.append( S_Iy[i] + h*Iy_1(t[i]+h/2, S_Ey[i]+(h/2)*Ey_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_Ey[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]) ) )
    S_Im.append( S_Im[i] + h*Im_1(t[i]+h/2, S_Em[i]+(h/2)*Em_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_Em[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]) ) )
    S_Io.append( S_Io[i] + h*Io_1(t[i]+h/2, S_Eo[i]+(h/2)*Eo_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_Eo[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]) ) )
    S_My.append( S_My[i] + h*My_1(t[i]+h/2, S_Iy[i]+(h/2)*Iy_1(t[i],S_Ey[i],S_Iy[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]) ) )
    S_Mm.append( S_Mm[i] + h*Mm_1(t[i]+h/2, S_Im[i]+(h/2)*Im_1(t[i],S_Em[i],S_Im[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]) ) )
    S_Mo.append( S_Mo[i] + h*Mo_1(t[i]+h/2, S_Io[i]+(h/2)*Io_1(t[i],S_Eo[i],S_Io[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]) ) )
    S_Dy.append( S_Dy[i] + h*Dy_1(t[i]+h/2, S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]) ) )
    S_Dm.append( S_Dm[i] + h*Dm_1(t[i]+h/2, S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]) ) )
    S_Do.append( S_Do[i] + h*Do_1(t[i]+h/2, S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]) ) )
    S_R.append( S_R[i] + h*R_1(t[i]+h/2, S_A[i]+(h/2)*A_1(t[i],S_Ey[i],S_Em[i],S_Eo[i],S_A[i]), S_My[i]+(h/2)*My_1(t[i],S_Iy[i],S_My[i]), S_Mm[i]+(h/2)*Mm_1(t[i],S_Im[i],S_Mm[i]), S_Mo[i]+(h/2)*Mo_1(t[i],S_Io[i],S_Mo[i]) ) )
    S_IT.append( S_IT[i] + h*IT_1(t[i]+h/2, S_Ey[i]+(h/2)*Ey_1(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy[i],S_Ey[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Em[i]+(h/2)*Em_1(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm[i],S_Em[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]), S_Eo[i]+(h/2)*Eo_1(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So[i],S_Eo[i],S_A[i],S_Iy[i],S_Im[i],S_Io[i],S_My[i],S_Mm[i],S_Mo[i],S_Dy[i],S_Dm[i],S_Do[i]) ) )


    S_Sy_x.append( S_Sy_x[i] + h*Sy_2(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sy_x[i]+(h/2)*Sy_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_Sm_x.append( S_Sm_x[i] + h*Sm_2(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sm_x[i]+(h/2)*Sm_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_So_x.append( S_So_x[i] + h*So_2(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_So_x[i]+(h/2)*So_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_Ey_x.append( S_Ey_x[i] + h*Ey_2(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sy_x[i]+(h/2)*Sy_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Ey_x[i]+(h/2)*Ey_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_Ey_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_Em_x.append( S_Em_x[i] + h*Em_2(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_Sm_x[i]+(h/2)*Sm_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Em_x[i]+(h/2)*Em_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_Em_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_Eo_x.append( S_Eo_x[i] + h*Eo_2(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]), S_So_x[i]+(h/2)*So_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Eo_x[i]+(h/2)*Eo_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_Eo_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]), S_Dy_x[i]+(h/2)*Dy_2(t[i],S_My_x[i]), S_Dm_x[i]+(h/2)*Dm_2(t[i],S_Mm_x[i]), S_Do_x[i]+(h/2)*Do_2(t[i],S_Mo_x[i]) ) )
    S_A_x.append( S_A_x[i] + h*A_2(t[i]+h/2, S_Ey_x[i]+(h/2)*Ey_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_Ey_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Em_x[i]+(h/2)*Em_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_Em_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Eo_x[i]+(h/2)*Eo_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_Eo_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]) ) )
    S_Iy_x.append( S_Iy_x[i] + h*Iy_2(t[i]+h/2, Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), S_Ey_x[i]+(h/2)*Ey_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_Ey_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]) ) )
    S_Im_x.append( S_Im_x[i] + h*Im_2(t[i]+h/2, Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), S_Em_x[i]+(h/2)*Em_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_Em_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]) ) )
    S_Io_x.append( S_Io_x[i] + h*Io_2(t[i]+h/2, Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), S_Eo_x[i]+(h/2)*Eo_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_Eo_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]) ) )
    S_My_x.append( S_My_x[i] + h*My_2(t[i]+h/2, Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), S_Iy_x[i]+(h/2)*Iy_2(t[i],Iy[i],S_Ey_x[i],S_Iy_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]) ) )
    S_Mm_x.append( S_Mm_x[i] + h*Mm_2(t[i]+h/2, Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), S_Im_x[i]+(h/2)*Im_2(t[i],Im[i],S_Em_x[i],S_Im_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]) ) )
    S_Mo_x.append( S_Mo_x[i] + h*Mo_2(t[i]+h/2, Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), S_Io_x[i]+(h/2)*Io_2(t[i],Io[i],S_Eo_x[i],S_Io_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]) ) )
    S_Dy_x.append( S_Dy_x[i] + h*Dy_2(t[i]+h/2, S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]) ) )
    S_Dm_x.append( S_Dm_x[i] + h*Dm_2(t[i]+h/2, S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]) ) )
    S_Do_x.append( S_Do_x[i] + h*Do_2(t[i]+h/2, S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]) ) )
    S_R_x.append( S_R_x[i] + h*R_2(t[i]+h/2, S_A_x[i]+(h/2)*A_2(t[i],S_Ey_x[i],S_Em_x[i],S_Eo_x[i],S_A_x[i]), S_My_x[i]+(h/2)*My_2(t[i],Iy[i],S_Iy_x[i],S_My_x[i]), S_Mm_x[i]+(h/2)*Mm_2(t[i],Im[i],S_Im_x[i],S_Mm_x[i]), S_Mo_x[i]+(h/2)*Mo_2(t[i],Io[i],S_Io_x[i],S_Mo_x[i]) ) )
    S_IT_x.append( S_IT_x[i] + h*IT_2(t[i]+h/2, S_Ey_x[i]+(h/2)*Ey_2(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sy_x[i],S_Ey_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Em_x[i]+(h/2)*Em_2(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_Sm_x[i],S_Em_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]), S_Eo_x[i]+(h/2)*Eo_2(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i],S_So_x[i],S_Eo_x[i],S_A_x[i],S_Iy_x[i],S_Im_x[i],S_Io_x[i],S_My_x[i],S_Mm_x[i],S_Mo_x[i],S_Dy_x[i],S_Dm_x[i],S_Do_x[i]) ) )






plt.plot(t,S_Sy_x,t,S_Sm_x,t,S_So_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of Susceptible Groups w.r.t parameter $x$")
plt.legend(["Sy","Sm","So"])
plt.grid(True)
plt.show()

plt.plot(t,S_Ey_x,t,S_Em_x,t,S_Eo_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of Exposed Groups w.r.t parameter $x$")
plt.legend(["Ey","Em","Eo"])
plt.grid(True)
plt.show()

plt.plot(t,S_Iy_x,t,S_Im_x,t,S_Io_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of Infectious Groups w.r.t parameter $x$")
plt.legend(["Iy","Im","Io"])
plt.grid(True)
plt.show()

plt.plot(t,S_My_x,t,S_Mm_x,t,S_Mo_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of Medical Treatment Groups w.r.t parameter $x$")
plt.legend(["My","Mm","Mo"])
plt.grid(True)
plt.show()

plt.plot(t,S_Dy_x,t,S_Dm_x,t,S_Do_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of Death Groups w.r.t parameter $x$")
plt.legend(["Dy","Dm","Do"])
plt.grid(True)
plt.show()

plt.plot(t,S_A_x,t,S_R_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities w.r.t parameter $x$")
plt.legend(["A","R"])
plt.grid(True)
plt.show()

plt.plot(t,S_IT,t,S_IT_x)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Sensitivities of IT(t) with logistics model")
plt.legend(["$\u03B2(t)$","x"])
plt.ylim(-450000,450000)
plt.grid(True)
plt.show()
