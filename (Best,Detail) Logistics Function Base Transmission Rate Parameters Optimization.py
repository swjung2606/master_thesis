from math import *
import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go
import plotly.io as pio
import requests
from lmfit import minimize, Parameters, Parameter, report_fit
from sympy import *


def f(xs,t,ps):
    try:
        a1 = ps["a1"].value
        a2 = ps["a2"].value
        a3 = ps["a3"].value
        k1 = ps["k1"].value
        k2 = ps["k2"].value
        b1 = ps["b1"].value
        b2 = ps["b2"].value
        epsilon = ps["epsilon"].value
        do = ps["do"].value
    except:
        a1,a2,a3,k1,k2,b1,b2,epsilon,do = ps

    Sy,Sm,So,Ey,Em,Eo,A,IT,ITy,ITm,ITo,Iy,Im,Io,My,Mm,Mo,Dy,Dm,Do = xs
    
    return [-( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) ,
            -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) ,
            -( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) ,
            ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sy*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Ey ,
            ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*Sm*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Em ,
            ( ((1.4081-b1)*a1)/((1.4081-b1)+(a1+b1-1.4081)*exp(k1*t)) + b1 + a3*a2/(a3+(a2-a3)*exp(-k2*t+b2)) )*(epsilon)*So*(Iy+Im+Io+0.01*A)*(1/(N-(My+Mm+Mo+Dy+Dm+Do))) - (eta)*Eo ,
            (1-p)*(eta)*(Ey+Em+Eo) - (gamma)*A , (eta)*p*(Ey+Em+Eo) , (eta)*p*Ey , (eta)*p*Em , (eta)*p*Eo , (eta)*p*Ey - x*Iy ,
            (eta)*p*Em - x*Im , (eta)*p*Eo - x*Io , x*Iy - ( alpha*(1-dy)+dy*theta )*My ,
            x*Im - ( alpha*(1-dm)+dm*theta )*Mm , x*Io - ( alpha*(1-do)+do*theta )*Mo , dy*(theta)*My , dm*(theta)*Mm , do*(theta)*Mo]


def g(t,ini0,ps):
    xx = odeint(f, ini0, t, args=(ps,))
    return xx


def residual(ps,ts,data,data_y,data_m,data_o,Death_m,Death_o):
    ini0 = 31451830,14578550,5194180, ps["init_Ey"].value,  ps["init_Em"].value, ps["init_Eo"].value, 0, 31,18,11,2,0,0,0,18,11,2,0,0,0
    model = g(ts,ini0,ps)
    return (model[:,7] - data).ravel() , (model[:,8] - data_y).ravel() , (model[:,9] - data_m).ravel() , (model[:,10] - data_o).ravel() , (model[:,18]-Death_m).ravel() , (model[:,19]-Death_o).ravel()


days=187
t=np.arange(0,days,1)
ini0 = [31451830,14578550,5194180,252,239,34,0,31,18,11,2,0,0,0,18,11,2,0,0,0]


N = 51225116
q = 0.01
eta = 1/5.1
p = 0.844
gamma = 1/14
alpha = 1/17.5
x = 1/4
dy = 0.0006
dm = 0.0104
do = 0.1344
theta = 1/13



a1=1.2694
a2=0.2631
a3=0.00011
k1=0.569
k2=0.05495
b1=0.14178
b2=-2.9
epsilon = 0.9953
    

data=[31,46,82,178,346,556,763,895,1147,1596,2023,2932,3527,4213,4813,5328,5766,6284,6767,7134,7382,7513,7755,7869,7979,8086,8162,8236,8320,8413,8565,8652,8799,8897,
          8961,9037,9137,9241,9332,9478,9583,9661,9786,9887,9976,10062,10156,10237,10284,10330,10384,10423,10450,10480,10512,10537,10564,10591,10612,10635,10653,10661,10674,10683,10694,10702,10708,10718,10728,10738,10752,10761,10765,10774,
          10780,10793,10801,10804,10806,10810,10822,10840,10874,10909,10936,10962,10991,11018,11037,11050,11065,11078,11110,11122,11142,11165,11190,11206,11225,11265,11344,11402,11441,11468,11503,11541,11590,11629,11668,11719,11776,11814,11852,
          11902,11947,12003,12051,12085,12121,12155,12198,12257,12306,12373,12421,12438,12484,12535,12563,12602,12653,12715,12757,12800,12851,12905,12968,13031,13092,13140,13184,13247,13297,13342,13377,13421,13483,13545,13584,13645,13705,13744,
          13778,13804,13849,13912,13971,14010,14123,14181,14206,14234,14282,14300,14336,14367,14397,14420,14454,14487,14530,14550,14593,14629,14657,14691,14745,14801,14904,15070,15349,15546,15792,16089,16377,16701,17033]

data_y=[18,27,44,84,155,270,386,481,630,931,1198,1770,2166,2576,2947,3263,3528,3796,4075,4286,4421,4492,4604,4662,4712,4758,4794,4826,4872,4911,4960,5000,5039,5087,5125,5166,5242
,5309,5352,5433,5506,5554,5607,5667,5724,5785,5837,5901,5933,5963,5995,6027,6042,6058,6080,6096,6114,6135,6149,6167,6180,6186,6196,6200,6209,6214,6217,6225,6235,6245,6255,6263,6267,6275
,6280,6288,6294,6296,6297,6301,6312,6329,6358,6392,6415,6437,6466,6488,6503,6514,6527,6533,6562,6574,6591,6601,6618,6627,6639,6661,6722,6760,6781,6801,6816,6828,6847,6859,6872,6895,6919
,6934,6943,6962,6987,7010,7024,7041,7064,7078,7105,7127,7149,7185,7206,7215,7240,7272,7287,7313,7346,7385,7409,7440,7461,7481,7519,7560,7594,7623,7651,7692,7727,7757,7774,7804,7853,7902
,7925,7970,8011,8042,8063,8085,8110,8143,8177,8193,8280,8312,8329,8348,8382,8394,8421,8438,8460,8473,8499,8520,8543,8556,8573,8593,8603,8618,8653,8689,8750,8840,8968,9058,9152,9283,9428,9577,9750]


data_m=[11,16,32,82,169,250,327,364,451,573,704,977,1140,1364,1549,1697,1826,1980,2117,2227,2307,2345,2432,2467,2508,2550,2580,2609,2635,2674,2722,2755,2790,2823,2841,2862,2878,2900,2922
,2981,3008,3030,3086,3110,3136,3153,3180,3193,3200,3212,3227,3231,3240,3253,3260,3267,3273,3276,3280,3284,3287,3288,3291,3294,3295,3296,3298,3300,3300,3300,3303,3304,3304,3304,3305,3307,3309
,3310,3311,3311,3312,3313,3317,3318,3321,3323,3323,3326,3330,3331,3333,3339,3341,3341,3343,3355,3362,3369,3374,3388,3406,3423,3437,3444,3462,3482,3506,3524,3541,3558,3581,3598,3617,3640,3654
,3684,3701,3716,3726,3740,3750,3772,3795,3820,3838,3845,3861,3876,3887,3897,3909,3928,3945,3954,3976,4002,4020,4037,4058,4072,4086,4101,4114,4125,4135,4147,4156,4165,4181,4193,4207,4214,4224
,4227,4235,4255,4276,4293,4315,4338,4341,4349,4358,4364,4371,4384,4390,4399,4405,4416,4433,4440,4463,4477,4494,4510,4524,4540,4574,4636,4750,4836,4964,5094,5196,5331,5453]


data_o=[2,3,6,12,22,36,50,50,66,92,121,185,221,273,317,368,412,508,575,621,654,676,719,740,759,778,788,801,813,828,883,897,970,987,995,1009,1017,1032,1058,1064,1069,1077,1093,1110,1116,1124
,1139,1143,1151,1155,1162,1165,1168,1169,1172,1174,1177,1180,1183,1184,1186,1187,1187,1189,1190,1192,1193,1193,1193,1193,1194,1194,1194,1195,1195,1198,1198,1198,1198,1198,1198,1198,1199,1199
,1200,1202,1202,1204,1204,1205,1205,1206,1207,1207,1208,1209,1210,1210,1212,1216,1216,1219,1223,1223,1225,1231,1237,1246,1255,1266,1276,1282,1292,1300,1306,1309,1326,1328,1331,1337,1343,1358,1362
,1368,1377,1378,1383,1387,1389,1392,1398,1402,1403,1406,1414,1422,1429,1434,1440,1445,1447,1454,1456,1460,1468,1470,1474,1478,1478,1482,1487,1488,1491,1492,1504,1514,1518,1524,1528,1531,1536
,1537,1542,1542,1544,1545,1547,1548,1550,1551,1554,1554,1557,1559,1560,1563,1568,1572,1580,1594,1631,1652,1676,1712,1753,1793,1830]



Death_y=[0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6]

Death_m=[0,0,0,0,1,1,3,3,4,4,6,8,8,9,10,12,13,16,16,16,17,18,20,20,20,20,20,20,22,22,24,24,24,24,25,28,30,30,
         31,31,31,31,32,33,34,36,37,38,39,40,40,41,43,45,46,47,47,47,47,48,49,
         49,49,49,50,50,50,50,50,50,50,50,50,50,50,51,51,51,52,52,52,52,52,52,
         52,53,53,53,53,53,53,53,53,53,54,54,54,54,54,54,54,54,54,54,54,54,54,
         54,54,54,54,54,54,54,55,55,55,55,55,55,55,55,55,55,55,56,56,56,56,56,
         56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57]


Death_o=[0,0,2,3,4,4,5,5,7,9,11,11,13,15,17,18,20,24,26,32,32,34,38,44,45,50,53,53,57,60,65,68,76,78,84,90,94,99,
         106,111,119,125,128,130,132,135,137,142,144,150,156,160,163,164,165,167,
         171,174,178,179,180,180,182,183,184,185,185,185,187,188,189,191,192,193,195,
         195,196,198,199,199,199,199,199,199,201,202,202,202,204,204,205,205,205,206,
         206,207,207,208,208,210,210,210,210,211,212,213,214,214,214,214,214,214,215,
         217,217,217,217,217,217,218,219,220,220,220,220,220,220,220,221,221,221,
         221,221,221,221,221,221,222,222,223,224,224,226,227,227,228,228,228,228,
         230,232,232,233,234,234,235,235,236,236,236,237,238,238,238,239,239,239,
         239,239,240,240,241,242,242,242,242,242,242,242,242,242,242,243,243,244,
         246,246]


params = Parameters()

params.add('init_Ey', value=float(ini0[3]), min=217.378, max=600.576)
params.add('init_Em', value=float(ini0[4]), min=203, max=600.2843)
params.add('init_Eo', value=float(ini0[5]), min=30, max=300.781)
params.add('a1', value=a1, min=1.1, max=1.3)
params.add('a2', value=a2, min=0.16, max=0.28)
params.add('a3', value=a3, min=0.00011, max=0.000113)
params.add('k1', value=k1, min=0.4, max=0.8)
params.add('k2', value=k2, min=0.048, max=0.0752)
params.add('b1', value=b1, min=0.1417, max=0.14185)
params.add('b2', value=b2, min=-3.9, max=-2)
params.add('epsilon', value=epsilon, min=0.701, max=1.7)
params.add('do', value=do, min=0.1, max=0.168)


result = minimize(residual, params, args=(t, data, data_y, data_m, data_o, Death_m,  Death_o), method='leastsq')
result.params
report_fit(result)
