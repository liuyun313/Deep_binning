
# coding: utf-8



# In[2]:

def usage ():
    print("usage:")
    print("  Python kmer.py -i inputfile -k k value")
    print("    -i inputfile should be in fasta format")
# In[1]:
import sys
import getopt

opts, args = getopt.getopt(sys.argv[1:], "hi:k:")
input_file=""
k=0
for op, value in opts:
    if op == "-i":
        input_file = value
    elif op == "-k":
        k=int(value)
    elif op == "-h":
        usage()
        sys.exit()


with open (input_file) as f:
    list=f.readlines()
    
N=0  #number of contigs 
### split contigs title and sequence into different lists
title=[]
contig=[]
a=""
title.append(list[0])
N+=1
for i in range(1,len(list)):
    list[i]=list[i].strip('\n')
    if(">" in list[i]):
        title.append(list[i])
        contig.append(a)
        a=""
        N+=1
    else:
       a+=list[i]

contig.append(a)        



# In[3]:

### convert ATGC to 0123 respectively
seq=[];
for i in range(0,len(contig)):
    x=str(contig[i]).replace("A","0").replace("T","1").replace("G","2").replace("C","3").replace(",","")
    seq.append(x)


# In[6]:


### compute k-mer frequency
import numpy as np

vec=np.zeros(pow(4,k)) # k^4，初始化vec向量
#mer=get_list(line,4) # k=4，获取k-mer片段存入mer 
weight=[]
for i in range(k):
    weight.append(pow(4,i)) # 4^(k-1),4^(k-2)……4^0，设定权重向量
weight=weight[::-1]

# 求取k-mer频率
nn=[]
#k=4
X=[]
for j in range(len(seq)):
    #mer[j]
    #mm=list(mer[j]) # 将一串字符串打散
    m=seq[j]
    for i in range(0,len(seq[j])-k+1): # 将字符串中的数字分别转换成整数并形成矩阵
        m1=m[i:i+k]
        a=0
        for ii in range(len(m1)):
            a+=weight[ii]*int(m1[ii])
        vec[a]=vec[a]+1 # 得到256维的k-mer频率向量vec
    # 反向互补
    for ii in range(len(vec)):
        if vec[ii]==-1:
            continue
        else:
            '''
            a=int(ii/64)
            b=int((ii-a*64)/16)
            c=int((ii-a*64-b*16)/4)
            d=ii-a*64-b*16-c*4 
            rank=np.array([a,b,c,d]) # 转化成4位二进制数
            '''
            ind=ii
            rank=np.zeros(k)
            for jj in range(k):
                rank[jj]=int(ind/pow(4,k-jj-1))
                ind-=rank[jj]*pow(4,k-jj-1)
            
            
            rank1=np.zeros(k)
            aa=np.argwhere(rank == 0) # 换成互补序列
            bb=np.argwhere(rank == 1)
            cc=np.argwhere(rank == 2)
            dd=np.argwhere(rank == 3)
            for jj in range(len(aa)):
                rank1[aa[jj]]=1
            for jj in range(len(bb)):
                rank1[bb[jj]]=0
            for jj in range(len(cc)):
                rank1[cc[jj]]=3
            for jj in range(len(dd)):
                rank1[dd[jj]]=2 
            rank1=rank1[::-1] # 倒序输出
            i1=0
            for jj in range(len(rank1)):
                i1+=rank1[jj]*weight[jj]
            if int(i1)!=ii:
                ee=vec[ii]+vec[int(i1)]
                vec[ii]=ee
                vec[int(i1)]=-1
    vecf=vec[vec!=-1] # 得到最后的k-mer频率矩阵vecf，共136维
    vec=np.zeros(pow(4,k))
    X.append(vecf)
#vecm=np.row_stack((vecm,vecf))
    
# In[]: normalization
    
def MaxMinNormalization(x,Max):
    x = x / Max;
    return x

X_=[]
for i in range(len(X)):
    x=MaxMinNormalization(X[i],np.max(X[i]))
    X_.append(x)
    


# In[19]:

import pandas as pd

#X1=np.array(X,dtype=float)
#kmer = [[float(x) for x in y.split(',') if len(x) >= 1] for y in X1[1:] if len(y) >= 1]
test=pd.DataFrame(index=title,data=X_)
k=str(k)
out_file=input_file+'_'+k+'mer'+'.csv'
test.to_csv(out_file,encoding='gbk',header=0)

'''
with open('kmer.csv', 'w',newline='') as csvfile:
    spamwriter = csv.writer(csvfile,dialect='excel')
    for i in range (len(title)):
        spamwriter.writerow([title[i],X1[i:]])
'''
print("Successful! Please go to kmer.csv for a look.")