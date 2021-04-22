import numpy as np
import random
from collections import Counter
import math
from matplotlib.pylab import plt
#delay=2#时延为2
#m=5#维度为，维度建议取值在3-7之间
def Permutation_Entropy(delay,m,data):
    def key(a):
        return a[0]
    X=[[]for i in range(len(data)-(m-1)*delay)]
    for i in range(len(X)):#将数据映射到m维空间
        for j in range(m):
            X[i].append([data[i+j*delay],j])
        X[i].sort(key=key)
    ordinal_patterns=[]
    for i in range(len(X)):
        s=''
        for j in range(m):
            s+=str(X[i][j][1])
        ordinal_patterns.append(s)

    ordinal_patterns=Counter(ordinal_patterns)
    P=[]
    for key in ordinal_patterns.keys():
        P.append(ordinal_patterns[key]/len(X))#计算序列模型概率分布
    H=0
    for i in range(len(P)):
        H+=P[i]*math.log(P[i])
    H*=-1.0
    return H/math.log(math.factorial(m))

import pandas as pd
df=pd.read_csv('data.csv')#加载数据
delay=1#时延为1
m=5#维度为5
print('时延为:',delay,'  维度为:',m)
for i in range(10):
    data=df.iloc[:,i]
    print('当前数据排列熵为：',Permutation_Entropy(1,5,data))
    plt.plot(data)
    plt.show()
