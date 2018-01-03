import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
''' 
##########################################################################
need complete last period 

t_min< number_samples/(total_time*frequency_object)<t_max
1.Provide t_min and t_max to run the program.
2.close the t_min and t_max to time period fastest the algorithm.
Obs-
3.Need higher sample rate
4.need samples with greater change during time each interval

Assuming total time =10 sec and total samples =1000
t_min< 100/frequency <t_max

########################################################################'''

"should be more variable with respect to time/change the limits on max & min  "
t_min=88
t_max=130
sample_rate =100
with open('05._0.9_1Hz.txt') as f:
    inp = []
    for line in f:
        line = line.split() # to deal with blank 
        if line:            # lines (ie skip them)
            line = [float(i) for i in line]
            inp.append(line)
#print inp;
for i in range(0,len(inp)):
        del inp[i][0]
inp = np.asarray(inp)   #converting the input list to array
#print inp
#plt.plot(inp)
#plt.show()
#print len(inp)

"---------------------------taking 2nd derivative ---------------------------"
def der(inp):
    d_plus1 = np.roll(inp,-1)
    d_plus2 = np.roll(inp,-2)
    d_plus3 = np.roll(inp,-3)
    d_minus1 = np.roll(inp,1)
    d_minus2 = np.roll(inp,2)
    d_minus3 = np.roll(inp,3)
    dd_inp =4*inp
    #print dd_inp
    dd_inp = (dd_inp+d_plus1+d_minus1-2*(d_plus2+d_minus2)-(d_plus3+d_minus3))/(16*t_s/samples*t_s/samples)
    #print dd_inp
    #plt.plot(dd_inp)
    #plt.show()
    return dd_inp

"-------------------Linear Warping by cubic spline interpolation--------------"
def warp(u,index):                         
    x =np.arange(len(u))
    tck = interpolate.splrep(x, u)    
    x1 =np.arange(index)
    return interpolate.splev(x1, tck)
"------------------segmented variance calculation ---------------------------"
def var_seg(seg,morph):
    warped = warp(morph,len(seg))
    var1 =np.subtract(seg,warped)
    return (np.linalg.norm(var1)*np.linalg.norm(var1))
   
"-----------------Update segmentation ---------------------------------------"
def segment(inp,u):
    s=[[0]*1]*(len(inp)+1)
    d=np.zeros(len(inp)+1)
    temp2=np.full(10*t_max-10*t_min+1,2e21)
    for i in range(1,len(inp)+1):
        temp=np.full(t_max-t_min+1,2e21)
        if t_min<=i<=t_max:
            d[i] = d[0]+var_seg(inp[0:i-1],u)

        if 2*t_min<=i<=2*t_max:
            if i-t_min<=t_max:
                for tau in range(t_min,i-t_min+1): 
                    #print temp,i,tau       #instead of this recursive fn check the time
                    temp[tau - t_min] = d[tau] + var_seg(inp[tau:i-1],u)
                d[i]=np.amin(temp)
                s[i]=np.append(s[i],int(np.argmin(temp)+t_min)) 
                #print i,s[i]
                temp=np.full(t_max-t_min+1,2e21)
            else:
                for tau in range(i-t_max,t_max+1):
                    #print temp,i,tau
                    temp[tau-i+t_max] = d[tau] + var_seg(inp[tau:i-1],u)
                d[i]=np.amin(temp)
                s[i]=np.append(s[i],int(np.argmin(temp)+i-t_max))
                #print i,s[i]
                temp=np.full(t_max-t_min+1,2e21)

        for itr in range(3,11):
                if itr*t_min<= i <=itr*t_max:
                    for tau in range(t_min,t_max+1):
                        if (itr-1)*t_min<=i-tau<=(itr-1)*t_max:
                            #print itr, temp,i,i-tau
                            temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
                    d[i]=np.amin(temp)
                    s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
                    print itr,i,s[i]
                    #print i, d[i]
                    temp=np.full(t_max-t_min+1,2e21)

                if i==len(inp) and itr*t_min< len(inp) < itr*t_max:
                    temp2 =np.full(10*t_max-10*t_min+1,2e21)
                    for last in range((itr-1)*t_max,len(inp)-1):
                        temp2[last-itr*t_min] = d[last] + var_seg(inp[last:i-1],u)
                        #print last, d[last], temp2[last-4*t_min]
                    d[i]=np.amin(temp2)
                    #print 'last',i,np.argmin(temp2)+4*t_min
                    s[len(inp)]=np.append(s[np.argmin(temp2)+itr*t_min],np.argmin(temp2)+itr*t_min)
                    #d[len(inp)]=np.amin(temp)
                    print itr
                    return s[i]
                    break

"----------------update template -------------------------------------------"
def template(inp,seg):
    prev_cost = 2e21
    cost = 0
    morph=[]
    for m in range(t_min,t_max+1):
        sums =np.zeros(m)
        for i in range(1,len(seg)):
            sums=sums+abs(seg[i]-seg[i-1])*warp(inp[seg[i-1]:seg[i]-1],m)
        sums=sums/len(inp)
        cost = cost_fn(inp,sums,seg)
        if cost<prev_cost:
            prev_cost =cost
            morph = sums
    return morph

"----------------variance for min cost-----------------------------"
def cost_fn(inp,u,s):
    sums=0
    for i in range(1,len(s)):
        sums=sums+var_seg(inp[s[i-1]:s[i]-1],u)
    return sums

"---------------iterations -------------------------------------------------"
def main(inp):
#    inp =der(inp)
    u = np.zeros(t_max)
    u=inp[0:t_max-1]
    s = []
    cost=np.zeros(15)
    prev_cost =2e21
    for i in range(0,1):
        #print i
        s =segment(inp,u)
        #print cost_fn(inp,u,s)
        print s 
        u =template(inp,s)
        #print u, len(u)
        cost[i+1]= cost_fn(inp,u,s)
        print cost[i+1],len(u)
        plt.plot(inp)
        #for t in range(0,len(s)):
        #    plt.axvline(s[t],color='r',linestyle='--')
        plt.plot(u)
        plt.show()
        if i>2 and (cost[i+1] == cost[i] or cost[i]==cost[i-2]):
            freq =0.000
            for j in range(1,len(s)):
                freq =freq + (s[j]-s[j-1])
            freq =freq/(len(s)-1)
            freq = sample_rate/freq
            break

        elif abs(cost[i+1] - cost[i]+1)/(cost[i]+cost[i+1]+1) < 0.5:
            freq =0.000
            for j in range(1,len(s)):
                freq =freq + (s[j]-s[j-1])
            freq =freq/(len(s)-1)
            #print freq
            freq = sample_rate/freq
            break
        #plt.figure(1+i) 
        #plt.plot(u)
    #plt.show()
    print freq,'Hz' 

main(inp)
