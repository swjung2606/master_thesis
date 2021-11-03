from math import *
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit



# South Korea 3rd model SEAIQR Simulation 


beta = 1.4081
beta1 = 0.4531
q = 0.01
eta = 1/5.1
p = 0.844
gamma = 1/14
alpha = 1/17.5

dy = 0.0006 # Based on 8/22 number of death
dm = 0.0104
do = 0.1312
theta = 1/13
epsilon = 1.1254
epsilon1 = 0.9957


# 2/18 - 2/29  : Beat -> Constant 1.4087

def Sy_prime(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -(beta)*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Sm_prime(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -(beta)*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def So_prime(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -(beta)*(epsilon1)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Ey_prime(t,Sy,Ey,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return (beta)*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey

def Em_prime(t,Sm,Em,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return (beta)*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em

def Eo_prime(t,So,Eo,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return (beta)*(epsilon1)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo

def A_prime(t,Ey,Em,Eo,A):
        return (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A

def AT_prime(t,Ey,Em,Eo):
        return (1-p)*(eta)*(Ey+Em+Eo)

def IT_prime(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def Iy_prime(t,Ey,Iy,x):
        return (eta)*p*Ey - (1/x)*Iy

def Im_prime(t,Em,Im,x):
        return (eta)*p*Em - (1/x)*Im

def Io_prime(t,Eo,Io,x):
        return (eta)*p*Eo - (1/x)*Io

def My_prime(t,Iy,My,x):
        return  (1/x)*Iy - ( alpha*(1-dy)+dy*theta )*My

def Mm_prime(t,Im,Mm,x):
        return  (1/x)*Im - ( alpha*(1-dm)+dm*theta )*Mm

def Mo_prime(t,Io,Mo,x):
        return  (1/x)*Io - ( alpha*(1-do)+do*theta )*Mo

def MT_prime(t,Iy,Im,Io,x):
        return (1/x)*(Iy+Im+Io)

def Dy_prime(t,My):
        return dy*(theta)*My

def Dm_prime(t,Mm):
        return dm*(theta)*Mm

def Do_prime(t,Mo):
        return do*(theta)*Mo

def R_prime(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )



# 3/1 - 8/22
def Sy_primea(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Sm_primea(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def So_primea(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*(epsilon1)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Ey_primea(t,Sy,Ey,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey

def Em_primea(t,Sm,Em,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em

def Eo_primea(t,So,Eo,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( (0.3771-0.0973)*sin((pi*(t-11))/350)+0.0973 )*(epsilon1)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo

def A_primea(t,Ey,Em,Eo,A):
        return (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A

def AT_primea(t,Ey,Em,Eo):
        return (1-p)*(eta)*(Ey+Em+Eo)

def IT_primea(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def Iy_primea(t,Ey,Iy,x):
        return (eta)*p*Ey - (1/x)*Iy

def Im_primea(t,Em,Im,x):
        return (eta)*p*Em - (1/x)*Im

def Io_primea(t,Eo,Io,x):
        return (eta)*p*Eo - (1/x)*Io

def My_primea(t,Iy,My,x):
        return  (1/x)*Iy - ( alpha*(1-dy)+dy*theta )*My

def Mm_primea(t,Im,Mm,x):
        return  (1/x)*Im - ( alpha*(1-dm)+dm*theta )*Mm

def Mo_primea(t,Io,Mo,x):
        return  (1/x)*Io - ( alpha*(1-do)+do*theta )*Mo

def MT_primea(t,Iy,Im,Io,x):
        return (1/x)*(Iy+Im+Io)

def Dy_primea(t,My):
        return dy*(theta)*My

def Dm_primea(t,Mm):
        return dm*(theta)*Mm

def Do_primea(t,Mo):
        return do*(theta)*Mo

def R_primea(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )




a1=1.2701
a2=0.1948
a3=0.00011
k1=0.7139
k2=0.0752
b1=0.1417
b2=-2.2807


# Logistics
def Sy_primel(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Sm_primel(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def So_primel(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Ey_primel(t,Sy,Ey,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey

def Em_primel(t,Sm,Em,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em

def Eo_primel(t,So,Eo,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo

def A_primel(t,Ey,Em,Eo,A):
        return (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A

def AT_primel(t,Ey,Em,Eo):
        return (1-p)*(eta)*(Ey+Em+Eo)

def IT_primel(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def Iy_primel(t,Ey,Iy,x):
        return (eta)*p*Ey - (1/x)*Iy

def Im_primel(t,Em,Im,x):
        return (eta)*p*Em - (1/x)*Im

def Io_primel(t,Eo,Io,x):
        return (eta)*p*Eo - (1/x)*Io

def My_primel(t,Iy,My,x):
        return  (1/x)*Iy - ( alpha*(1-dy)+dy*theta )*My

def Mm_primel(t,Im,Mm,x):
        return  (1/x)*Im - ( alpha*(1-dm)+dm*theta )*Mm

def Mo_primel(t,Io,Mo,x):
        return  (1/x)*Io - ( alpha*(1-do)+do*theta )*Mo

def MT_primel(t,Iy,Im,Io):
        return (1/x)*(Iy+Im+Io)

def Dy_primel(t,My):
        return dy*(theta)*My

def Dm_primel(t,Mm):
        return dm*(theta)*Mm

def Do_primel(t,Mo):
        return do*(theta)*Mo

def R_primel(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )






t=np.linspace(0,186,187)
h=t[1]
dd=[]
ee=[]
x=np.linspace(2,5,31)

# 2/18 - 2/29

for j in range(0,31):
        Sy=[31451830]
        Sm=[14578550]
        So=[5194180]
        Ey=[252]
        Em=[239]
        Eo=[34]
        A=[0]
        AT=[0]
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
        R0=[]
        IT=[31]
        MT=[31]

        Syl=[31451865]
        Sml=[14578212]
        Sol=[5194184]
        Eyl=[217]
        Eml=[577]
        Eol=[30]
        Al=[0]
        ATl=[0]
        Iyl=[0]
        Iml=[0]
        Iol=[0]
        Myl=[18]
        Mml=[11]
        Mol=[2]
        Dyl=[0]
        Dml=[0]
        Dol=[0]
        Rl=[0]
        R0l=[]
        ITl=[31]
        MTl=[31]




        N = 51225116
        for i in range(0,11):
                Sy.append( Sy[i] + h*Sy_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                Sm.append( Sm[i] + h*Sm_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                So.append( So[i] + h*So_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                Ey.append( Ey[i] + h*Ey_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                Em.append( Em[i] + h*Em_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                Eo.append( Eo[i] + h*Eo_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
                A.append( A[i] + h*A_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]) ) )
                AT.append( AT[i] + h*AT_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
                IT.append( IT[i] + h*IT_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
                Iy.append( Iy[i] + h*Iy_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), x[j] ) )
                Im.append( Im[i] + h*Im_prime(t[i]+h/2, Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), x[j] ) )
                Io.append( Io[i] + h*Io_prime(t[i]+h/2, Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), x[j] ) )
                My.append( My[i] + h*My_prime(t[i]+h/2, Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i],x[j]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), x[j] ) )
                Mm.append( Mm[i] + h*Mm_prime(t[i]+h/2, Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), x[j] ) )
                Mo.append( Mo[i] + h*Mo_prime(t[i]+h/2, Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]), x[j] ) )
                Dy.append( Dy[i] + h*Dy_prime(t[i]+h/2, My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]) ) )
                Dm.append( Dm[i] + h*Dm_prime(t[i]+h/2, Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]) ) )
                Do.append( Do[i] + h*Do_prime(t[i]+h/2, Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]) ) )   
                R.append( R[i] + h*R_prime(t[i]+h/2, A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i],x[j]) ) )




        for i in range(11,186):
                Sy.append( Sy[i] + h*Sy_primea(t[i]+h/2, Sy[i]+(h/2)*Sy_primea(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                Sm.append( Sm[i] + h*Sm_primea(t[i]+h/2, Sm[i]+(h/2)*Sm_primea(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                So.append( So[i] + h*So_primea(t[i]+h/2, So[i]+(h/2)*So_primea(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                Ey.append( Ey[i] + h*Ey_primea(t[i]+h/2, Sy[i]+(h/2)*Sy_primea(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Ey[i]+(h/2)*Ey_primea(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                Em.append( Em[i] + h*Em_primea(t[i]+h/2, Sm[i]+(h/2)*Sm_primea(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primea(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                Eo.append( Eo[i] + h*Eo_primea(t[i]+h/2, So[i]+(h/2)*So_primea(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primea(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), Dy[i]+(h/2)*Dy_primea(t[i],My[i]), Dm[i]+(h/2)*Dm_primea(t[i],Mm[i]), Do[i]+(h/2)*Do_primea(t[i],Mo[i]) ) )
                A.append( A[i] + h*A_primea(t[i]+h/2, Ey[i]+(h/2)*Ey_primea(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primea(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primea(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]) ) )
                AT.append( AT[i] + h*AT_primea(t[i]+h/2, Ey[i]+(h/2)*Ey_primea(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primea(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primea(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
                IT.append( IT[i] + h*IT_primea(t[i]+h/2, Ey[i]+(h/2)*Ey_primea(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primea(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primea(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
                Iy.append( Iy[i] + h*Iy_primea(t[i]+h/2, Ey[i]+(h/2)*Ey_primea(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), x[j] ) )
                Im.append( Im[i] + h*Im_primea(t[i]+h/2, Em[i]+(h/2)*Em_primea(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), x[j] ) )
                Io.append( Io[i] + h*Io_primea(t[i]+h/2, Eo[i]+(h/2)*Eo_primea(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), x[j] ) )
                My.append( My[i] + h*My_primea(t[i]+h/2, Iy[i]+(h/2)*Iy_primea(t[i],Ey[i],Iy[i],x[j]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), x[j] ) )
                Mm.append( Mm[i] + h*Mm_primea(t[i]+h/2, Im[i]+(h/2)*Im_primea(t[i],Em[i],Im[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), x[j] ) )
                Mo.append( Mo[i] + h*Mo_primea(t[i]+h/2, Io[i]+(h/2)*Io_primea(t[i],Eo[i],Io[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]), x[j] ) )
                Dy.append( Dy[i] + h*Dy_primea(t[i]+h/2, My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]) ) )
                Dm.append( Dm[i] + h*Dm_primea(t[i]+h/2, Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]) ) )
                Do.append( Do[i] + h*Do_primea(t[i]+h/2, Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]) ) )   
                R.append( R[i] + h*R_primea(t[i]+h/2, A[i]+(h/2)*A_primea(t[i],Ey[i],Em[i],Eo[i],A[i]), My[i]+(h/2)*My_primea(t[i],Iy[i],My[i],x[j]), Mm[i]+(h/2)*Mm_primea(t[i],Im[i],Mm[i],x[j]), Mo[i]+(h/2)*Mo_primea(t[i],Io[i],Mo[i],x[j]) ) )

        for i in range(0,186):
                Syl.append( Syl[i] + h*Sy_primel(t[i]+h/2, Syl[i]+(h/2)*Sy_primel(t[i],Syl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Sml.append( Sml[i] + h*Sm_primel(t[i]+h/2, Sml[i]+(h/2)*Sm_primel(t[i],Sml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Sol.append( Sol[i] + h*So_primel(t[i]+h/2, Sol[i]+(h/2)*So_primel(t[i],Sol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Eyl.append( Eyl[i] + h*Ey_primel(t[i]+h/2, Syl[i]+(h/2)*Sy_primel(t[i],Syl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eyl[i]+(h/2)*Ey_primel(t[i],Syl[i],Eyl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Eml.append( Eml[i] + h*Em_primel(t[i]+h/2, Sml[i]+(h/2)*Sm_primel(t[i],Sml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eml[i]+(h/2)*Em_primel(t[i],Sml[i],Eml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Eol.append( Eol[i] + h*Eo_primel(t[i]+h/2, Sol[i]+(h/2)*So_primel(t[i],Sol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eol[i]+(h/2)*Eo_primel(t[i],Sol[i],Eol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), Dyl[i]+(h/2)*Dy_primel(t[i],Myl[i]), Dml[i]+(h/2)*Dm_primel(t[i],Mml[i]), Dol[i]+(h/2)*Do_primel(t[i],Mol[i]) ) )
                Al.append( Al[i] + h*A_primel(t[i]+h/2, Eyl[i]+(h/2)*Ey_primel(t[i],Syl[i],Eyl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eml[i]+(h/2)*Em_primel(t[i],Sml[i],Eml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eol[i]+(h/2)*Eo_primel(t[i],Sol[i],Eol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]) ) )
                ATl.append( ATl[i] + h*AT_primel(t[i]+h/2, Eyl[i]+(h/2)*Ey_primel(t[i],Syl[i],Eyl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eml[i]+(h/2)*Em_primel(t[i],Sml[i],Eml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eol[i]+(h/2)*Eo_primel(t[i],Sol[i],Eol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]) ) )
                ITl.append( ITl[i] + h*IT_primel(t[i]+h/2, Eyl[i]+(h/2)*Ey_primel(t[i],Syl[i],Eyl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eml[i]+(h/2)*Em_primel(t[i],Sml[i],Eml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Eol[i]+(h/2)*Eo_primel(t[i],Sol[i],Eol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]) ) )
                Iyl.append( Iyl[i] + h*Iy_primel(t[i]+h/2, Eyl[i]+(h/2)*Ey_primel(t[i],Syl[i],Eyl[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), x[j] ) )
                Iml.append( Iml[i] + h*Im_primel(t[i]+h/2, Eml[i]+(h/2)*Em_primel(t[i],Sml[i],Eml[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), x[j] ) )
                Iol.append( Iol[i] + h*Io_primel(t[i]+h/2, Eol[i]+(h/2)*Eo_primel(t[i],Sol[i],Eol[i],Al[i],Iyl[i],Iml[i],Iol[i],Myl[i],Mml[i],Mol[i],Dyl[i],Dml[i],Dol[i]), Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), x[j] ) )
                Myl.append( Myl[i] + h*My_primel(t[i]+h/2, Iyl[i]+(h/2)*Iy_primel(t[i],Eyl[i],Iyl[i],x[j]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), x[j] ) )
                Mml.append( Mml[i] + h*Mm_primel(t[i]+h/2, Iml[i]+(h/2)*Im_primel(t[i],Eml[i],Iml[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), x[j] ) )
                Mol.append( Mol[i] + h*Mo_primel(t[i]+h/2, Iol[i]+(h/2)*Io_primel(t[i],Eol[i],Iol[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]), x[j] ) )
                Dyl.append( Dyl[i] + h*Dy_primel(t[i]+h/2, Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]) ) )
                Dml.append( Dml[i] + h*Dm_primel(t[i]+h/2, Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]) ) )
                Dol.append( Dol[i] + h*Do_primel(t[i]+h/2, Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]) ) )   
                Rl.append( Rl[i] + h*R_primel(t[i]+h/2, Al[i]+(h/2)*A_primel(t[i],Eyl[i],Eml[i],Eol[i],Al[i]), Myl[i]+(h/2)*My_primel(t[i],Iyl[i],Myl[i],x[j]), Mml[i]+(h/2)*Mm_primel(t[i],Iml[i],Mml[i],x[j]), Mol[i]+(h/2)*Mo_primel(t[i],Iol[i],Mol[i],x[j]) ) )


        dd.append(IT[186])
        ee.append(ITl[186])
                
     




data=[31,46,82,178,346,556,763,895,1147,1596,2023,2932,3527,4213,4813,5328,5766,6284,6767,7134,7382,7513,7755,7869,7979,8086,8162,8236,8320,8413,8565,8652,8799,8897,
          8961,9037,9137,9241,9332,9478,9583,9661,9786,9887,9976,10062,10156,10237,10284,10330,10384,10423,10450,10480,10512,10537,10564,10591,10612,10635,10653,10661,10674,10683,10694,10702,10708,10718,10728,10738,10752,10761,10765,10774,
          10780,10793,10801,10804,10806,10810,10822,10840,10874,10909,10936,10962,10991,11018,11037,11050,11065,11078,11110,11122,11142,11165,11190,11206,11225,11265,11344,11402,11441,11468,11503,11541,11590,11629,11668,11719,11776,11814,11852,
          11902,11947,12003,12051,12085,12121,12155,12198,12257,12306,12373,12421,12438,12484,12535,12563,12602,12653,12715,12757,12800,12851,12905,12968,13031,13092,13140,13184,13247,13297,13342,13377,13421,13483,13545,13584,13645,13705,13744,
          13778,13804,13849,13912,13971,14010,14123,14181,14206,14234,14282,14300,14336,14367,14397,14420,14454,14487,14530,14550,14593,14629,14657,14691,14745,14801,14904,15070,15349,15546,15792,16089,16377,16701,17033]
         
Death_y=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6]

Death_m=[12,13,16,16,16,17,18,20,20,20,20,20,20,22,22,24,24,24,24,25,28,30,30,
         31,31,31,31,32,33,34,36,37,38,39,40,40,41,43,45,46,47,47,47,47,48,49,
         49,49,49,50,50,50,50,50,50,50,50,50,50,50,51,51,51,52,52,52,52,52,52,
         52,53,53,53,53,53,53,53,53,53,54,54,54,54,54,54,54,54,54,54,54,54,54,
         54,54,54,54,54,54,54,55,55,55,55,55,55,55,55,55,55,55,56,56,56,56,56,
         56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57]

Death_o=[18,20,24,26,32,32,34,38,44,45,50,53,53,57,60,65,68,76,78,84,90,94,99,
         106,111,119,125,128,130,132,135,137,142,144,150,156,160,163,164,165,167,
         171,174,178,179,180,180,182,183,184,185,185,185,187,188,189,191,192,193,195,
         195,196,198,199,199,199,199,199,199,201,202,202,202,204,204,205,205,205,206,
         206,207,207,208,208,210,210,210,210,211,212,213,214,214,214,214,214,214,215,
         217,217,217,217,217,217,218,219,220,220,220,220,220,220,220,221,221,221,
         221,221,221,221,221,221,222,222,223,224,224,226,227,227,228,228,228,228,
         230,232,232,233,234,234,235,235,236,236,236,237,238,238,238,239,239,239,
         239,239,240,240,241,242,242,242,242,242,242,242,242,242,242,243,243,244,
         246,246]
       



def func(x,a,b):
	return a*np.exp(b*x)

popt1,pcov1 = curve_fit(func,x,dd)
popt2,pcov2 = curve_fit(func,x,ee)




plt.plot(x,func(x,*popt1),color="blue")
plt.plot(x,func(x,*popt2),color="orange")
plt.scatter(x,dd,marker="o",color="blue")
plt.scatter(x,ee,marker="*",color="orange")
plt.xlabel("parameter $1/x$")
plt.title("Simulation of total cumulative cases on Aug 22")
plt.legend(["Sine Function-Based Model Fitting","Logistics Funcgion-Based Model Fitting","Total cumuative cases on Aug 22 (Sine model)","Total cumulative cases on Aug 22(Logistics model)"])
plt.text(2,90000,"Fitting Curve = 0.4777 * exp(2.5613(1/x))",color="Blue")
plt.text(2,60000,"Fitting Curve = 0.4995 * exp(2.5494(1/x))",color="orange")
plt.grid(True)

plt.show()

