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
    tr = x[1]
    a = x[2]
    b = x[3]
    c = (1 - 1 / b)
    return np.sum((tr + (ts - tr) / ((1 + (a * w) ** b) ** c) - t) ** 2)


# The constraint of the coefficient should be greater than zero
def constraint1(x):
    ts = x[0]
    tr = x[1]
    a = x[2]
    b = x[3]
    c = (1 - 1 / b)
    f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
    beta = (len(t) * np.sum(t * f) - np.sum(t) * np.sum(f)) / (len(f) * np.sum(t ** 2) - np.sum(t) ** 2)
    return beta - 1.


def constraint2(x):
    ts = x[0]
    tr = x[1]
    a = x[2]
    b = x[3]
    c = (1 - 1 / b)
    f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
    Rsquared = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
    return Rsquared - 1.


# The constraint and bound
con1 = {'type': 'eq', 'fun': constraint1}
con2 = {'type': 'eq', 'fun': constraint2}
cons1 = [con1, con2]
cons2 = [con1]
cons3 = [con2]
bond1 = (0, 1)
bond2 = (0, 1)
bond3 = (0, 100)
bond4 = (0, 1000)
bond = [bond1, bond2, bond3, bond4]
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
TR=[]
A=[]
B=[]
RS=[]

print("The optimization started please wait")
fragment=1000
i=0
while i<fragment:
    TS_r=random.randint(int(t[0]*100), 100)/100
    TR_r=random.randint(0, int(t[len(t)-1]*100+0.2*(t[0]-t[len(t)-1])*100))/100
    A_r=(100 * random.randint(1, 9) +10 * random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10 + random.randint(1, 9) / 100) / (
                    10 ** random.randint(0, 6))
    B_r=(100 * random.randint(1, 9) +10 * random.randint(1, 9) + random.randint(1, 9) + random.randint(1, 9) / 10 + random.randint(1, 9) / 100) / (
                    10 ** random.randint(1, 3))
    try:
        x = [TS_r, TR_r, A_r, B_r]
        ts = x[0]
        tr = x[1]
        a = x[2]
        b = x[3]
        c = (1 - 1 / b)
        f_limit=tr + (ts - tr) / ((1 + (a * 10**6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        if r >= 0.5 and r <= 1:
            TS.append(ts)
            TR.append(tr)
            A.append(a)
            B.append(b)
            RS.append(r)
            i = i + 1
        else:i=i

    except:
        do_nothing=0

maxiteration=20
TSF=[]
TRF=[]
AF=[]
BF=[]
RSF=[]
for i in range(0,len(TS)):
    x = [TS[i], TR[i], A[i], B[i]]
    sol4 = minimize(objective, x, method='SLSQP', bounds=bond,
                    options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
    try:
        ts = sol4.x[0]
        tr = sol4.x[1]
        a = sol4.x[2]
        b = sol4.x[3]
        c = (1 - 1 / b)
        f_limit = tr + (ts - tr) / ((1 + (a * 10 ** 6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        TRF.append(tr)
        AF.append(a)
        BF.append(b)
        RSF.append(r)

    except:
        do_nothing=0


J=0
while J==0:
    print("Elapsed time=" + str(time.process_time() - start))
    old_r = np.mean(RS)
    new_r = np.mean(RSF)
    print("previous R^2=" + str(old_r))
    print("current R^2=" + str(new_r))
    if old_r<new_r and new_r-old_r>0.00001:
        TS = TSF
        TR = TRF
        A = AF
        B = BF
        RS = RSF

        TSF = []
        TRF = []
        AF = []
        BF = []
        RSF = []
        for zz in range(0, len(TS)):
            x = [TS[zz], TR[zz], A[zz], B[zz]]
            sol4 = minimize(objective, x, method='SLSQP', bounds=bond,
                            options={'maxiter': maxiteration, 'ftol': 1e-100000000000000, 'disp': False})
            try:
                ts = sol4.x[0]
                tr = sol4.x[1]
                a = sol4.x[2]
                b = sol4.x[3]
                c = (1 - 1 / b)
                f_limit = tr + (ts - tr) / ((1 + (a * 10 ** 6) ** b) ** c)
                f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
                r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
                if r>=new_r:
                    TSF.append(ts)
                    TRF.append(tr)
                    AF.append(a)
                    BF.append(b)
                    RSF.append(r)
                else:do_nothing=0


            except:
                do_nothing = 1
    else: J=4


"""clustring"""
if len(TSF)>1:
    selected_fragment_array = np.stack((TSF, TRF, AF, BF), axis=-1)
    i = 1
    J = 0
    while i <= len(TSF) and J == 0:
        RSQCL = []
        kmeans = KMeans(n_clusters=i, random_state=0).fit(selected_fragment_array)
        y = kmeans.cluster_centers_
        i = i + 1
        TS = y[:, 0]
        TR = y[:, 1]
        A = y[:, 2]
        B = y[:, 3]
        for zz in range(0, len(TS)):
            try:
                ts = TS[zz]
                tr = TR[zz]
                a = A[zz]
                b = B[zz]
                c = (1 - 1 / b)
                f_limit = tr + (ts - tr) / ((1 + (a * 10 ** 6) ** b) ** c)
                f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
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
maxiteration=1000
TSF=[]
TRF=[]
AF=[]
BF=[]
ERF=[]
RSF=[]
for i in range(0,len(TS)):
    x = [TS[i], TR[i], A[i], B[i]]
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
        tr = sol1.x[1]
        a = sol1.x[2]
        b = sol1.x[3]
        c = (1 - 1 / b)
        f_limit=tr + (ts - tr) / ((1 + (a * 10**6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        TRF.append(tr)
        AF.append(a)
        BF.append(b)
        ERF.append(np.sum((tr + (ts - tr) / ((1 + (a * w) ** b) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol2.x[0]
        tr = sol2.x[1]
        a = sol2.x[2]
        b = sol2.x[3]
        c = (1 - 1 / b)
        f_limit=tr + (ts - tr) / ((1 + (a * 10**6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        TRF.append(tr)
        AF.append(a)
        BF.append(b)
        ERF.append(np.sum((tr + (ts - tr) / ((1 + (a * w) ** b) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol3.x[0]
        tr = sol3.x[1]
        a = sol3.x[2]
        b = sol3.x[3]
        c = (1 - 1 / b)
        f_limit=tr + (ts - tr) / ((1 + (a * 10**6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        TRF.append(tr)
        AF.append(a)
        BF.append(b)
        ERF.append(np.sum((tr + (ts - tr) / ((1 + (a * w) ** b) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i
    try:
        ts = sol4.x[0]
        tr = sol4.x[1]
        a = sol4.x[2]
        b = sol4.x[3]
        c = (1 - 1 / b)
        f_limit=tr + (ts - tr) / ((1 + (a * 10**6) ** b) ** c)
        f = tr + (ts - tr) / ((1 + (a * w) ** b) ** c)
        r = 1 - (np.sum((t - f) ** 2)) / (np.sum((t - np.mean(t)) ** 2))
        TSF.append(ts)
        TRF.append(tr)
        AF.append(a)
        BF.append(b)
        ERF.append(np.sum((tr + (ts - tr) / ((1 + (a * w) ** b) ** c) - t) ** 2))
        RSF.append(r)

    except:
        i=i

min_index = int(ERF.index(min(ERF)))
print("-------------------------results---------------------------")
print("Theta s="+str(TSF[min_index]))
print("Theta r="+str(TRF[min_index]))
print("a="+str(AF[min_index]))
print("n="+str(BF[min_index]))
print("sum of squared error="+str(ERF[min_index]))
print("R^2="+str(RSF[min_index]))
