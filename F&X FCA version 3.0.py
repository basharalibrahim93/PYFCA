#the DATA file should be located in the program directory and label as "DATA.dat"
#first column should contain the data label
#second column should contain the suction pressure in kPa
#third column should contain the dimensionless volumetric water content (not in percentage)
#only tabs should exist between the first, second, and third columns( no comma)

import numpy as np
from scipy.optimize import minimize
np.seterr(all=None, divide='ignore', over='ignore', under=None, invalid='ignore')
import time
import random
from sklearn.cluster import KMeans
start = time.process_time()


# The objective function which is used to minimize the sum of square error
def objective(x):
    ts = x[0]
    a = x[1]
    b = x[2]
    c = x[3]
    wr = x[4]
    return np.sum(((1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
            (np.log(np.exp(1) + (w / a) ** b)) ** c) - t) ** 2)


# The constraint of the coefficient should be greater than zero
def constraint1(x):
    ts = x[0]
    a = x[1]
    b = x[2]
    c = x[3]
    wr = x[4]
    f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
    beta = (len(t) * np.sum(t * f) - np.sum(t) * np.sum(f)) / (len(f) * np.sum(t ** 2) - np.sum(t) ** 2)
    return beta - 1.


def constraint2(x):
    ts = x[0]
    a = x[1]
    b = x[2]
    c = x[3]
    wr = x[4]
    f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
    Rsquared = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
    return Rsquared - 1.


# The constraint and bound
con1 = {'type': 'eq', 'fun': constraint1}
con2 = {'type': 'eq', 'fun': constraint2}
cons1 = [con1, con2]
cons2 = [con1]
cons3 = [con2]
bond1 = (0, 1)
bond2 = (0, 10**5)
bond3 = (0, 1000)
bond4 = (0, 1000)
bond5 = (0, 1000)
bond5 = (0, 10**6)
bond = [bond1, bond2, bond3, bond4,bond5 ]
Frgmnt = []
Sol = []
SolOPT = []
inputfile = np.genfromtxt("DATA.dat")
data_num = np.array(inputfile[0:, 0])
data_w = np.array(inputfile[0:, 1])
data_t = np.array(inputfile[0:, 2])
w = data_w
t = data_t
FN=[]
TS=[]
A=[]
B=[]
C=[]
WR=[]
RS=[]

print("The optimization started please wait")
fragment=1000
i=0
while i<fragment:
    TS_r=random.randint(int(t[0]*100), 100)/100
    A_r=(1000 * random.randint(1, 9) +100 * random.randint(1, 9) + 10*random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10) / (
                    10 ** random.randint(0, 11))
    B_r=(100 * random.randint(1, 9) +10 * random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10 + random.randint(1, 9) / 100) / (
                    10 ** random.randint(0, 16))
    C_r=(100 * random.randint(1, 9) +10 * random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10 + random.randint(1, 9) / 100) / (
                    10 ** random.randint(0, 4))
    W_r=(100 * random.randint(1, 9) +10 * random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10 + random.randint(1, 9) / 100) *(
                    10 ** random.randint(-1, 4))
    if W_r>=10**6:W_r=99999
    try:
        x = [TS_r, A_r, B_r,C_r,W_r]
        ts = x[0]
        a = x[1]
        b = x[2]
        c = x[3]
        wr = x[4]
        f_limit= (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (10**6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        if r >= 0.5 and r <= 1:
            TS.append(ts)
            WR.append(wr)
            A.append(a)
            B.append(b)
            C.append(c)
            RS.append(r)
            i = i + 1
        else:i=i

    except:
        do_nothing=0

maxiteration=20
TSF=[]
WRF=[]
AF=[]
BF=[]
CF=[]
RSF=[]
for i in range(0,len(TS)):
    x = [TS[i],  A[i], B[i], C[i],WR[i]]
    sol4 = minimize(objective, x, method='SLSQP', bounds=bond,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    try:
        ts = sol4.x[0]
        a = sol4.x[1]
        b = sol4.x[2]
        c = sol4.x[3]
        wr = sol4.x[4]
        f_limit= (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (10**6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        WRF.append(wr)
        AF.append(a)
        BF.append(b)
        CF.append(c)
        RSF.append(r)

    except:
        do_nothing=0


J=0
while J==0:
    print("Elapsed time="+str(time.process_time() - start))
    old_r=np.mean(RS)
    new_r=np.mean(RSF)
    print("previous R^2="+str(old_r))
    print("current R^2="+str(new_r))
    if old_r<new_r and new_r-old_r>0.00001:
        TS = TSF
        WR = WRF
        A = AF
        B = BF
        C = CF
        RS = RSF

        TSF = []
        WRF = []
        AF = []
        BF = []
        CF = []
        RSF = []
        for zz in range(0, len(TS)):
            x = [TS[zz], A[zz], B[zz], C[zz], WR[zz]]
            sol4 = minimize(objective, x, method='SLSQP', bounds=bond,
                            options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
            try:
                ts = sol4.x[0]
                a = sol4.x[1]
                b = sol4.x[2]
                c = sol4.x[3]
                wr = sol4.x[4]
                f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                            (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
                f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
                r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
                if r>=new_r:
                    TSF.append(ts)
                    WRF.append(wr)
                    AF.append(a)
                    BF.append(b)
                    CF.append(c)
                    RSF.append(r)
                else:do_nothing=0


            except:
                do_nothing = 1
    else: J=4


"""clustring"""
if len(TSF)>1:
    selected_fragment_array = np.stack((TSF, AF, BF, CF, WRF), axis=-1)
    i = 1
    J = 0
    while i <= len(TSF) and J == 0:
        RSQCL = []
        kmeans = KMeans(n_clusters=i, random_state=0).fit(selected_fragment_array)
        y = kmeans.cluster_centers_
        i = i + 1
        TS = y[:, 0]
        A = y[:, 1]
        B = y[:, 2]
        C = y[:, 3]
        WR = y[:, 4]
        for zz in range(0, len(TS)):
            try:
                ts = TS[zz]
                wr = WR[zz]
                a = A[zz]
                b = B[zz]
                c = C[zz]
                f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                        (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
                f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
                r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
                RSQCL.append(r)
            except:
                RSQCL.append(0)
        print("number of cluster: " + str(i - 1) + " R^2= " + str(np.mean(RSQCL)))
        if np.mean(RSQCL) >= np.mean(RSF):
            J = 2


        else:
            J = 0


else:print("number of cluster:0 ")

maxiteration=200
TSF=[]
WRF=[]
AF=[]
BF=[]
CF=[]
ERF=[]
RSF=[]
for i in range(0,len(TS)):
    x = [TS[i], A[i], B[i],C[i], WR[i]]
    sol1 = minimize(objective, x, method='SLSQP', bounds=bond, constraints=cons1,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    sol2 = minimize(objective, x, method='SLSQP', bounds=bond, constraints=cons2,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    sol3 = minimize(objective, x, method='SLSQP', bounds=bond, constraints=cons3,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    sol4 = minimize(objective, x, method='SLSQP', bounds=bond,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    try:
        ts = sol1.x[0]
        a = sol1.x[1]
        b = sol1.x[2]
        c = sol1.x[3]
        wr = sol1.x[4]
        f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        WRF.append(wr)
        AF.append(a)
        BF.append(b)
        CF.append(c)
        ERF.append(np.sum(((1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
            (np.log(np.exp(1) + (w / a) ** b)) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol2.x[0]
        a = sol2.x[1]
        b = sol2.x[2]
        c = sol2.x[3]
        wr = sol2.x[4]
        f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        WRF.append(wr)
        AF.append(a)
        BF.append(b)
        CF.append(c)
        ERF.append(np.sum(((1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
            (np.log(np.exp(1) + (w / a) ** b)) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol3.x[0]
        a = sol3.x[1]
        b = sol3.x[2]
        c = sol3.x[3]
        wr = sol3.x[4]
        f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        WRF.append(wr)
        AF.append(a)
        BF.append(b)
        CF.append(c)
        ERF.append(np.sum(((1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
            (np.log(np.exp(1) + (w / a) ** b)) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol4.x[0]
        a = sol4.x[1]
        b = sol4.x[2]
        c = sol4.x[3]
        wr = sol4.x[4]
        f_limit = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
                (np.log(np.exp(1) + (10 ** 6 / a) ** b)) ** c)
        f = (1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / ((np.log(np.exp(1) + (w / a) ** b)) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        WRF.append(wr)
        AF.append(a)
        BF.append(b)
        CF.append(c)
        ERF.append(np.sum(((1 - np.log(1 + w / wr) / np.log(1 + 10 ** 6 / wr)) * ts / (
            (np.log(np.exp(1) + (w / a) ** b)) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
min_index = int(ERF.index(min(ERF)))
print("-------------------------results---------------------------")
print("Theta s="+str(TSF[min_index]))
print("a="+str(AF[min_index]))
print("n="+str(BF[min_index]))
print("m="+str(CF[min_index]))
print("residual suction="+str(WRF[min_index]))
print("sum of squared error="+str(ERF[min_index]))
print("R^2="+str(RSF[min_index]))