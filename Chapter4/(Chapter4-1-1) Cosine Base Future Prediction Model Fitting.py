from math import *
import numpy as np
import matplotlib.pyplot as plt
import csv




# South Korea 3rd model SEAIQR Simulation 

Sy=[31451865]
Sm=[14578212]
So=[5194184]
Ey=[217]
Em=[577]
Eo=[30]
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
bbt=[]


IT=[31]
ITy=[18]
ITm=[11]
ITo=[2]
aa=[]


N = 51225116
beta = 1.4081
beta1 = 0.4531
q = 0.01
eta = 1/5.1
p = 0.844
gamma = 1/14
alpha = 1/17.5
x = 1/4 # 1/4.1 looks best...but I computed it 1/4
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


# 2/18 - 8/22

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

def AT_prime(t,Ey,Em,Eo):
        return (1-p)*(eta)*(Ey+Em+Eo)

def IT_prime(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def ITy_prime(t,Ey):
        return (eta)*p*Ey

def ITm_prime(t,Em):
        return (eta)*p*Em

def ITo_prime(t,Eo):
        return (eta)*p*Eo

def Iy_prime(t,Ey,Iy):
        return (eta)*p*Ey - x*Iy

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

def MT_prime(t,Iy,Im,Io):
        return x*(Iy+Im+Io)

def Dy_prime(t,My):
        return dy*(theta)*My

def Dm_prime(t,Mm):
        return dm*(theta)*Mm

def Do_prime(t,Mo):
        return do*(theta)*Mo

def R_prime(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )




# 8/23 -

c1=0.3847
c2=0.2961
l1=60

def Sy_primef(t,Sy,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Sm_primef(t,Sm,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def So_primef(t,So,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return -( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) 

def Ey_primef(t,Sy,Ey,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*Sy*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey

def Em_primef(t,Sm,Em,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*Sm*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em

def Eo_primef(t,So,Eo,A,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do):
        return ( 0.5*(c1-c2)*cos( (t-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) )*(epsilon)*So*(Iy+Im+Io+q*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo

def A_primef(t,Ey,Em,Eo,A):
        return (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A

def AT_primef(t,Ey,Em,Eo):
        return (1-p)*(eta)*(Ey+Em+Eo)

def IT_primef(t,Ey,Em,Eo):
        return (eta)*p*(Ey+Em+Eo)

def ITy_primef(t,Ey):
        return (eta)*p*Ey

def ITm_primef(t,Em):
        return (eta)*p*Em

def ITo_primef(t,Eo):
        return (eta)*p*Eo

def Iy_primef(t,Ey,Iy):
        return (eta)*p*Ey - x*Iy

def Im_primef(t,Em,Im):
        return (eta)*p*Em - x*Im

def Io_primef(t,Eo,Io):
        return (eta)*p*Eo - x*Io

def My_primef(t,Iy,My):
        return  x*Iy - ( alpha*(1-dy)+dy*theta )*My

def Mm_primef(t,Im,Mm):
        return  x*Im - ( alpha*(1-dm)+dm*theta )*Mm

def Mo_primef(t,Io,Mo):
        return  x*Io - ( alpha*(1-do)+do*theta )*Mo

def Dy_primef(t,My):
        return dy*(theta)*My

def Dm_primef(t,Mm):
        return dm*(theta)*Mm

def Do_primef(t,Mo):
        return do*(theta)*Mo

def R_primef(t,A,My,Mm,Mo):
        return (gamma)*A + (alpha)*( (1-dy)*My + (1-dm)*Mm + (1-do)*Mo )




t=np.linspace(0,478,479)
t1=np.linspace(0,317,318)
h=t[1]




for i in range(0,187):  
    Sy.append( Sy[i] + h*Sy_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Sm.append( Sm[i] + h*Sm_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    So.append( So[i] + h*So_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Ey.append( Ey[i] + h*Ey_prime(t[i]+h/2, Sy[i]+(h/2)*Sy_prime(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Em.append( Em[i] + h*Em_prime(t[i]+h/2, Sm[i]+(h/2)*Sm_prime(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    Eo.append( Eo[i] + h*Eo_prime(t[i]+h/2, So[i]+(h/2)*So_prime(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_prime(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_prime(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_prime(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_prime(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_prime(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_prime(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_prime(t[i],My[i]), Dm[i]+(h/2)*Dm_prime(t[i],Mm[i]), Do[i]+(h/2)*Do_prime(t[i],Mo[i]) ) )
    A.append( A[i] + h*A_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_prime(t[i],Ey[i],Em[i],Eo[i],A[i]) ) )
    AT.append( AT[i] + h*AT_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
    IT.append( IT[i] + h*IT_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
    ITy.append( ITy[i] + h*ITy_prime(t[i]+h/2, Ey[i]+(h/2)*Ey_prime(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
    ITm.append( ITm[i] + h*ITm_prime(t[i]+h/2, Em[i]+(h/2)*Em_prime(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
    ITo.append( ITo[i] + h*ITo_prime(t[i]+h/2, Eo[i]+(h/2)*Eo_prime(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
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


for i in range(187,478):
        
        Sy.append( Sy[i] + h*Sy_primef(t[i]+h/2, Sy[i]+(h/2)*Sy_primef(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        Sm.append( Sm[i] + h*Sm_primef(t[i]+h/2, Sm[i]+(h/2)*Sm_primef(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        So.append( So[i] + h*So_primef(t[i]+h/2, So[i]+(h/2)*So_primef(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        Ey.append( Ey[i] + h*Ey_primef(t[i]+h/2, Sy[i]+(h/2)*Sy_primef(t[i],Sy[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        Em.append( Em[i] + h*Em_primef(t[i]+h/2, Sm[i]+(h/2)*Sm_primef(t[i],Sm[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        Eo.append( Eo[i] + h*Eo_primef(t[i]+h/2, So[i]+(h/2)*So_primef(t[i],So[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]), Dy[i]+(h/2)*Dy_primef(t[i],My[i]), Dm[i]+(h/2)*Dm_primef(t[i],Mm[i]), Do[i]+(h/2)*Do_primef(t[i],Mo[i]) ) )
        A.append( A[i] + h*A_primef(t[i]+h/2, Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]) ) )
        AT.append( AT[i] + h*AT_primef(t[i]+h/2, Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
        IT.append( IT[i] + h*IT_primef(t[i]+h/2, Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
        ITy.append( ITy[i] + h*ITy_primef(t[i]+h/2, Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
        ITm.append( ITm[i] + h*ITm_primef(t[i]+h/2, Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
        ITo.append( ITo[i] + h*ITo_primef(t[i]+h/2, Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]) ) )
        Iy.append( Iy[i] + h*Iy_primef(t[i]+h/2, Ey[i]+(h/2)*Ey_primef(t[i],Sy[i],Ey[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]) ) )
        Im.append( Im[i] + h*Im_primef(t[i]+h/2, Em[i]+(h/2)*Em_primef(t[i],Sm[i],Em[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]) ) )
        Io.append( Io[i] + h*Io_primef(t[i]+h/2, Eo[i]+(h/2)*Eo_primef(t[i],So[i],Eo[i],A[i],Iy[i],Im[i],Io[i],My[i],Mm[i],Mo[i],Dy[i],Dm[i],Do[i]), Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]) ) )
        My.append( My[i] + h*My_primef(t[i]+h/2, Iy[i]+(h/2)*Iy_primef(t[i],Ey[i],Iy[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]) ) )
        Mm.append( Mm[i] + h*Mm_primef(t[i]+h/2, Im[i]+(h/2)*Im_primef(t[i],Em[i],Im[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]) ) )
        Mo.append( Mo[i] + h*Mo_primef(t[i]+h/2, Io[i]+(h/2)*Io_primef(t[i],Eo[i],Io[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]) ) )
        Dy.append( Dy[i] + h*Dy_primef(t[i]+h/2, My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]) ) )
        Dm.append( Dm[i] + h*Dm_primef(t[i]+h/2, Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]) ) )
        Do.append( Do[i] + h*Do_primef(t[i]+h/2, Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]) ) )   
        R.append( R[i] + h*R_primef(t[i]+h/2, A[i]+(h/2)*A_primef(t[i],Ey[i],Em[i],Eo[i],A[i]), My[i]+(h/2)*My_primef(t[i],Iy[i],My[i]), Mm[i]+(h/2)*Mm_primef(t[i],Im[i],Mm[i]), Mo[i]+(h/2)*Mo_primef(t[i],Io[i],Mo[i]) ) )




for i in range(0,187):
    aa.append( ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t[i])) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t[i]+b2)) ) )

for i in range(187,479):
    aa.append( ( 0.5*(c1-c2)*cos( (t[i]-187+(l1/pi)*acos((0.673-c1-c2)/(c1-c2)))*(pi/l1) ) + 0.5*(c1+c2) ) )





data=[31,46,82,178,346,556,763,895,1147,1596,2023,2932,3527,4213,4813,5328,5766,6284,6767,7134,7382,7513,7755,7869,7979,8086,8162,8236,8320,8413,8565,8652,8799,8897,
          8961,9037,9137,9241,9332,9478,9583,9661,9786,9887,9976,10062,10156,10237,10284,10330,10384,10423,10450,10480,10512,10537,10564,10591,10612,10635,10653,10661,10674,10683,10694,10702,10708,10718,10728,10738,10752,10761,10765,10774,
          10780,10793,10801,10804,10806,10810,10822,10840,10874,10909,10936,10962,10991,11018,11037,11050,11065,11078,11110,11122,11142,11165,11190,11206,11225,11265,11344,11402,11441,11468,11503,11541,11590,11629,11668,11719,11776,11814,11852,
          11902,11947,12003,12051,12085,12121,12155,12198,12257,12306,12373,12421,12438,12484,12535,12563,12602,12653,12715,12757,12800,12851,12905,12968,13031,13092,13140,13184,13247,13297,13342,13377,13421,13483,13545,13584,13645,13705,13744,
          13778,13804,13849,13912,13971,14010,14123,14181,14206,14234,14282,14300,14336,14367,14397,14420,14454,14487,14530,14550,14593,14629,14657,14691,14745,14801,14904,15070,15349,15546,15792,16089,16377,16701,17033,
          17430,17696,17976,18296,18737,19108,19431,19730,19978,20213,20480,20675,20873,21041,21208,21327,21463,21619,21774,21950,22086,22207,22316,22422,22535,22688,22814,22924,23006,23076,23137,23247,23372,23486,23547,23642,23692,23730,
          23843,23920,23983,24058,24122,24195,24270,24384,24453,24507,24579,24637,24735,24837,24921,25031,25078,25149,25238,25322,25377,25467,25586,25740,25810,25866,25979,
          26065,26164,26278,26384,26498,26617,26707,26779,26897,27022,27169,27258,27401,27572,27627,27773,27916,28107,28312,28520,28743,28973,29286,29629,29992,30378,30708,30979,31328,31710,
          32293,32862,33366,33816,34254,34705,35216,35756,36835,36968,37599,38214,38808,39494,40176,40865,41815,42845,43563,44443,45521,46599,47661,48714,
49811,50737,51606,52698,53683,54924,56056,57026,57834,58880,59930,60897]





data_y=[18,27,44,84,155,270,386,481,630,931,1198,1770,2166,2576,2947,3263,3528,3796,4075,4286,4421,4492,4604,4662,4712,4758,4794,4826,4872,4911,4960,5000,5039,5087,5125,5166,5242
,5309,5352,5433,5506,5554,5607,5667,5724,5785,5837,5901,5933,5963,5995,6027,6042,6058,6080,6096,6114,6135,6149,6167,6180,6186,6196,6200,6209,6214,6217,6225,6235,6245,6255,6263,6267,6275
,6280,6288,6294,6296,6297,6301,6312,6329,6358,6392,6415,6437,6466,6488,6503,6514,6527,6533,6562,6574,6591,6601,6618,6627,6639,6661,6722,6760,6781,6801,6816,6828,6847,6859,6872,6895,6919
,6934,6943,6962,6987,7010,7024,7041,7064,7078,7105,7127,7149,7185,7206,7215,7240,7272,7287,7313,7346,7385,7409,7440,7461,7481,7519,7560,7594,7623,7651,7692,7727,7757,7774,7804,7853,7902
,7925,7970,8011,8042,8063,8085,8110,8143,8177,8193,8280,8312,8329,8348,8382,8394,8421,8438,8460,8473,8499,8520,8543,8556,8573,8593,8603,8618,8653,8689,8750,8840,8968,9058,9152,9283,9428,9577,9750,
9937,10063,10201,10353,10561,10728,10876,11018,11125,11229,11343,11429,11513,11574,11644,11698,11757,11810,11868,11937,11992,12048,12092,12142,12186,12255,12316,12344,12380,12419,12448,12489
,12540,12589,12615,12662,12683,12709,12767,12807,12838,12870,12900,12937,12986,13020,13053,13086,13114,13151,13215,13284,13347,13393,13421,13450,13494,13542,13575,13621,13671,13732,13764,13800
,13873,13915,13973,14031,14091,14151,14233,14283,14325,14385,14458,14548,14600,14664,14732,14790,14851,14941,15028,15151,15252,15363,15485,15661,15877,16096,16360,16573,16748,16962,17192,17571
,17902,18199,18481,18723,18996,19354,19691,20033,20352,20639,20917,21202,21531,21905,22251,22742,23221,23598,24029,24542,25055,25550,26050,26600,27029,27449,27927,28429,29067,29648
,30146,30542,31093,31580,32055]



data_m=[11,16,32,82,169,250,327,364,451,573,704,977,1140,1364,1549,1697,1826,1980,2117,2227,2307,2345,2432,2467,2508,2550,2580,2609,2635,2674,2722,2755,2790,2823,2841,2862,2878,2900,2922
,2981,3008,3030,3086,3110,3136,3153,3180,3193,3200,3212,3227,3231,3240,3253,3260,3267,3273,3276,3280,3284,3287,3288,3291,3294,3295,3296,3298,3300,3300,3300,3303,3304,3304,3304,3305,3307,3309
,3310,3311,3311,3312,3313,3317,3318,3321,3323,3323,3326,3330,3331,3333,3339,3341,3341,3343,3355,3362,3369,3374,3388,3406,3423,3437,3444,3462,3482,3506,3524,3541,3558,3581,3598,3617,3640,3654
,3684,3701,3716,3726,3740,3750,3772,3795,3820,3838,3845,3861,3876,3887,3897,3909,3928,3945,3954,3976,4002,4020,4037,4058,4072,4086,4101,4114,4125,4135,4147,4156,4165,4181,4193,4207,4214,4224
,4227,4235,4255,4276,4293,4315,4338,4341,4349,4358,4364,4371,4384,4390,4399,4405,4416,4433,4440,4463,4477,4494,4510,4524,4540,4574,4636,4750,4836,4964,5094,5196,5331,5453,5622,5722,5825
,5942,6115,6251,6372,6500,6604,6695,6794,6864,6941,7021,7095,7138,7198,7273,7347,7421,7484,7533,7583,7626,7673,7736,7784,7845,7874,7896,7922,7964,8020,8061,8087,8126,8148,8157,8201,8227,8251
,8278,8298,8324,8344,8403,8426,8441,8474,8486,8510,8535,8552,8574,8589,8625,8661,8681,8699,8728,8778,8830,8858,8873,8906,8942,8971,9013,9049,9083,9111,9139,9165,9208,9252,9295,9317,9363,9397
,9425,9476,9516,9585,9653,9742,9835,9917,10013,10112,10217,10310,10398,10468,10571,10689,10856,11041,11202,11331,11477,11614,11716,11886,12098,12281,12515,12753,12952,13203,13428,13653,13975
,14350,14596,14880,15286,15692,16103,16497,16918,17274,17586,18021,18375,18828,19208,19578,19879,20259,20661,21021]



data_o=[2,3,6,12,22,36,50,50,66,92,121,185,221,273,317,368,412,508,575,621,654,676,719,740,759,778,788,801,813,828,883,897,970,987,995,1009,1017,1032,1058,1064,1069,1077,1093,1110,1116,1124
,1139,1143,1151,1155,1162,1165,1168,1169,1172,1174,1177,1180,1183,1184,1186,1187,1187,1189,1190,1192,1193,1193,1193,1193,1194,1194,1194,1195,1195,1198,1198,1198,1198,1198,1198,1198,1199,1199
,1200,1202,1202,1204,1204,1205,1205,1206,1207,1207,1208,1209,1210,1210,1212,1216,1216,1219,1223,1223,1225,1231,1237,1246,1255,1266,1276,1282,1292,1300,1306,1309,1326,1328,1331,1337,1343,1358,1362
,1368,1377,1378,1383,1387,1389,1392,1398,1402,1403,1406,1414,1422,1429,1434,1440,1445,1447,1454,1456,1460,1468,1470,1474,1478,1478,1482,1487,1488,1491,1492,1504,1514,1518,1524,1528,1531,1536
,1537,1542,1542,1544,1545,1547,1548,1550,1551,1554,1554,1557,1559,1560,1563,1568,1572,1580,1594,1631,1652,1676,1712,1753,1793,1830,1871,1911,1950,2001,2061,2129,2183,2212,2249,2289,2343,2382
,2419,2446,2469,2491,2508,2536,2559,2592,2610,2626,2641,2654,2676,2697,2714,2735,2752,2761,2767,2794,2812,2836,2845,2854,2861,2864,2875,2886,2894,2910,2924,2934,2940,2961,2974,2980,2991,3000
,3010,3018,3022,3064,3068,3074,3083,3099,3103,3118,3137,3178,3188,3193,3200,3208,3220,3234,3244,3264,3273,3285,3289,3304,3312,3326,3341,3374,3398,3412,3446,3459,3494,3508,3526,3545,3571,3612
,3640,3679,3708,3737,3763,3795,3829,3866,3919,3965,4004,4054,4095,4146,4179,4254,4335,4445,4544,4654,4760,4843,4961,5098,5274,5369,5534,5693,5852,6008,6167,6293,6434,6571,6750,6879
,7029,7200,7302,7413,7528,7689,7821]

 
Death_y=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
         6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,9
,9,9,9,10,10,10,10,10,10,10]



Death_m=[12,13,16,16,16,17,18,20,20,20,20,20,20,22,22,24,24,24,24,25,28,30,30,
         31,31,31,31,32,33,34,36,37,38,39,40,40,41,43,45,46,47,47,47,47,48,49,
         49,49,49,50,50,50,50,50,50,50,50,50,50,50,51,51,51,52,52,52,52,52,52,
         52,53,53,53,53,53,53,53,53,53,54,54,54,54,54,54,54,54,54,54,54,54,54,
         54,54,54,54,54,54,54,55,55,55,55,55,55,55,55,55,55,55,56,56,56,56,56,
         56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,58,58,58,58,58,58,
         58,58,58,58,58,58,58,58,59,59,59,59,59,59,59,59,59,59,60,60,61,61,61,
         62,63,63,63,63,63,63,64,65,65,65,66,66,66,66,67,68,68,68,68,68,68,68,
         70,70,71,73,74,74,74,74,74,74,74,74,74,74,76,76,76,76,76,76,76,77,78,
         78,78,80,80,80,80,80,82,83,84,84,85,85,86,86,86,87,87,88,88,88,88,88,
         88,88,88,89,90,90,93,93,95,97,99,99,103,103,105,109,113,118,119,119
,120,122,124,126,129,133,136]



Death_o=[18,20,24,26,32,32,34,38,44,45,50,53,53,57,60,65,68,76,78,84,90,94,99,
         106,111,119,125,128,130,132,135,137,142,144,150,156,160,163,164,165,167,
         171,174,178,179,180,180,182,183,184,185,185,185,187,188,189,191,192,193,195,
         195,196,198,199,199,199,199,199,199,201,202,202,202,204,204,205,205,205,206,
         206,207,207,208,208,210,210,210,210,211,212,213,214,214,214,214,214,214,215,
         217,217,217,217,217,217,218,219,220,220,220,220,220,220,220,221,221,221,
         221,221,221,221,221,221,222,222,223,224,224,226,227,227,228,228,228,228,
         230,232,232,233,234,234,235,235,236,236,236,237,238,238,238,239,239,239,
         239,239,240,240,241,242,242,242,242,242,242,242,242,242,242,243,243,244,
         246,246,246,246,247,249,250,253,257,259,260,260,262,265,267,269,270,272,
         277,280,282,286,290,293,298,302,302,307,312,313,318,320,322,322,326,328,
         332,333,337,338,344,346,347,351,351,351,351,354,355,356,358,360,360,360,
         364,365,367,369,370,370,371,374,376,377,377,377,380,381,382,383,384,386,
         388,390,392,393,394,395,396,398,402,403,403,404,406,407,408,408,410,410,
         412,413,415,418,419,421,423,424,429,430,432,432,432,434,441,445,450,450
,454,456,459,467,472,478,478,483,494,505,523,534,546,557,576,595,611,628,643
,661,674,683,720,736,754]




# RK4 method IT
IT_RK4 = [31, 145.13590548292382, 260.4530622682586, 393.7886014598257, 558.1737842483321, 765.4202521928521, 1026.623518495136, 1350.5679059363497, 1739.921391350226, 2186.777799620483, 2671.782482646318, 3169.785385204217, 3658.361745194712, 4122.876152967507, 4556.523526042067, 4957.8840530672405, 5328.402833239129, 5670.679991774325, 5987.543941485101, 6281.636701693977, 6555.275445304383, 6810.442446746873, 7048.824530921201, 7271.864661106718, 7480.810219683786, 7676.753021379577, 7860.660618627468, 8033.400106370437, 8195.75600873624, 8348.443705399444, 8492.119580288643, 8627.388793508982, 8754.811338811518, 8874.90686344386, 8988.158589437588, 9095.01657579972, 9195.900490297865, 9291.20200975809, 9381.286933023595, 9466.497066518423, 9547.151925536582, 9623.550282677777, 9695.971586697737, 9764.677269336466, 9829.911953672416, 9891.90457470075, 9950.86942078943, 10007.007103179749, 10060.50545959813, 10111.540397216471, 10160.276679558681, 10206.86866144598, 10251.460975664959, 10294.189174703872, 10335.180330616009, 10374.55359582147, 10412.420727441347, 10448.886577564943, 10484.049551676666, 10518.00203731145, 10550.830804863492, 10582.61738234105, 10613.438405738534, 10643.365946585158, 10672.46781812599, 10700.807861495694, 10728.446213156853, 10755.439554793002, 10781.841346770801, 10807.702046215662, 10833.069310680148, 10857.988188324223, 10882.501295470372, 10906.64898234462, 10930.46948776586, 10953.99908350064, 10977.272208958133, 11000.32159686022, 11023.178390484163, 11045.872253040136, 11068.431469712434, 11090.883042861704, 11113.25278085551, 11135.565380966142, 11157.84450674751, 11180.112860277226, 11202.392249625533, 11224.7036518894, 11247.067272108066, 11269.502598355246, 11292.02845328336, 11314.66304237629, 11337.423999149392, 11360.328427518665, 11383.392941545286, 11406.633702746765, 11430.066455152195, 11453.706558265972, 11477.569018092312, 11501.668516361497, 11526.019438088339, 11550.635897583525, 11575.531763029485, 11600.72067972401, 11626.216092087127, 11652.03126451953, 11678.179301194315, 11704.673164857617, 11731.525694708163, 11758.749623420612, 11786.357593372733, 11814.362172132174, 11842.775867254499, 11871.611140440495, 11900.880421097356, 11930.596119345202, 11960.770638507554, 11991.416387121699, 12022.545790502472, 12054.171301890761, 12086.305413215916, 12118.960665499426, 12152.149658925406, 12185.885062601861, 12220.179624035207, 12255.046178339117, 12290.497657197528, 12326.547097600453, 12363.207650370145, 12400.492588494175, 12438.41531528103, 12476.989372352993, 12516.228447490239, 12556.146382339375, 12596.757179998904, 12638.075012493511, 12680.114228148424, 12722.889358874563, 12766.415127374676, 12810.706454280175, 12855.778465227932, 12901.646497885868, 12948.326108935802, 12995.833081021632, 13044.183429670593, 13093.393410195036, 13143.479524581842, 13194.458528376324, 13246.34743756723, 13299.163535479192, 13352.924379678743, 13407.64780889984, 13463.351949994605, 13520.055224914842, 13577.776357729686, 13636.534381684613, 13696.348646306873, 13757.238824562295, 13819.224920068247, 13882.327274367462, 13946.5665742673, 14011.963859248934, 14078.540528950844, 14146.318350730935, 14215.319467311494, 14285.56640451114, 14357.08207906789, 14429.889806557323, 14504.013309409875, 14579.47672503116, 14656.304614029232, 14734.521968552634, 14814.154220743045, 14895.227251306347, 14977.767398205851, 15061.801465481458, 15147.356732198477, 15234.46096152983, 15323.14240997536, 15413.42983672195, 15505.352513148171, 15598.940232477167, 15694.223319581499, 15791.232640943668, 15889.999614776068, 15990.556221304107, 16092.935013216267]
       

plt.subplot(1,2,1)
plt.plot(t,Iy,t,Im,t,Io,marker="*")
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.legend(["$I_{y}$(0~49)","$I_{m}$(50~69)","$I_{o}$(70~)"])
plt.title("Number of Infectious Compartments")
plt.grid(True)


plt.subplot(1,2,2)
plt.scatter(t1,data,s=5,color="red")
plt.plot(t,IT,color="blue")
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Total Cumulative Cases")
plt.legend(["Simulation","Observed Data up to Dec 31, 2020"])
plt.grid(True)


plt.show()


####################################


plt.scatter(t1,data_y, s=2, color="blue")
plt.scatter(t1,data_m, s=2, color="orange")
plt.scatter(t1,data_o, s=2, color="green")
plt.plot(t,ITy,t,ITm,t,ITo)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.title("Cumulative Cases by age groups")
plt.legend(["Culumative Cases Simulation (Young)","Culumative Cases Simulation (Middle)","Culumative Cases Simulation (Old)","Observed Culumative Cases (Young)","Observed Culumative Cases (Middle)","Observed Culumative Cases (Old)"])
plt.grid(True)
plt.show()


####################################


plt.scatter(t1[15:],Death_y, s=2, color="blue")
plt.scatter(t1[15:],Death_m, s=2, color="orange")
plt.scatter(t1[15:],Death_o, s=2, color="green")
plt.plot(t,Dy,t,Dm,t,Do)
plt.xlabel("Days (t, t=0 is Feb 18)")
plt.ylabel("Cumulative Death Cases")
plt.legend(["$D_{y}$(0~49)","$D_{m}$(50~69)","$D_{o}$(70~)","Observed Death(Y)","Observed Death(M)","Observed Death(O)"])
plt.grid(True)

plt.show()





####################################


fig, ax_left = plt.subplots()
ax_right = ax_left.twinx()

ax_left.plot(t,Iy,t,Im,t,Io,marker="*")
ax_right.plot(t,aa, color='red')
plt.title("Active Cases and Transmission Rate")
plt.legend(["Transmission Rate"])
plt.grid(True)
plt.show()
