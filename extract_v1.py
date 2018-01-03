import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
''' 
##########################################################################
not used 'for loop' for iterations of segmentation 

******This is working fine with incomplete period data, but need favourable data

t_min< number_samples/(total_time*frequency_object)<t_max
1.Provide t_min and t_max to run the program.
2.close the t_min and t_max to time period fastest the algorithm.
Obs-
3.Need higher sample rate
4.need samples with greater change during time each interval

Assuming total time =10 sec and total samples =1000
t_min< 100/frequency <t_max

#########################################################################
'''
t_min=130
t_max=156
sample_rate =100
with open('0.5_0.7_1Hz.txt') as f:
    inp = []
    for line in f:
        line = line.split() # to deal with blank 
        if line:            # lines (ie skip them)
            line = [float(i) for i in line]
            inp.append(line)
#print inp;
for i in range(0,len(inp)):
        del inp[i][1]
inp = np.asarray(inp)   #converting the input list to array
#min_value= np.amin(inp)
#inp = np.subtract(inp,min_value)
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
    #plt.plot(t_range,dd_inp)
    #plt.show()
    return dd_inp

"-------------------Linear Warping by cubic spline interpolation--------------"
def warp(u,index):                              #have to add max and min. limit on index
    x =np.arange(len(u))
    #tck = interpolate.splrep(x, u,seg_min,seg_max)
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
    
    #print s
    d=np.zeros(len(inp)+1)
    for i in range(1,len(inp)+1):
        temp=np.full(t_max-t_min+1,1e20)
        if t_min<=i<=t_max:
            d[i] = d[0]+var_seg(inp[0:i-1],u)

        if 2*t_min<=i<=2*t_max:
            if i-t_min<=t_max:
                for tau in range(t_min,i-t_min+1): 
                    #print temp,i,tau       #instead of this recursive fn check the time
                    temp[tau - t_min] = d[tau] + var_seg(inp[tau:i-1],u)
                d[i]=np.amin(temp)
                s[i]=np.append(s[i],int(np.argmin(temp)+t_min))
                #print np.argsort(temp)+t_min, i,'2' 
                #print i,s[i] ,d[i]
                temp=np.full(t_max-t_min+1,1e20)
            #if not(i-t_max>=t_min and i-t_min==t_max):
            else:
                for tau in range(i-t_max,t_max+1):
                    #print temp,i,tau
                    temp[tau-i+t_max] = d[tau] + var_seg(inp[tau:i-1],u)
                d[i]=np.amin(temp)
                s[i]=np.append(s[i],int(np.argmin(temp)+i-t_max))
                #print np.argsort(temp)+t_min, i,'2'
                #print i,s[i], d[i]
                temp=np.full(t_max-t_min+1,1e20)

        elif 3*t_min<= i <=3*t_max:
            for tau in range(t_min,t_max+1):
                if 2*t_min<=i-tau<=2*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            #print i-np.argsort(temp)-t_min, i,'3'
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i], d[i]
            temp=np.full(t_max-t_min+1,1e20)

        elif 4*t_min<= i<= 4*t_max:
            for tau in range(t_min,t_max+1):
                if 3*t_min<=i-tau<=3*t_max:
                    #print temp,i,i-tau, '4th temp cal' 
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            #print d[i],'4th'
            #print i-np.argsort(temp)-t_min, i,'4'
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i], d[i],'4th'
            temp=np.full(t_max-t_min+1,1e20)

        elif 5*t_min<=i<=5*t_max:
            for tau in range(t_min,t_max+1):
                if 4*t_min<=i-tau<=4*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i], d[i],'5th'
            temp=np.full(t_max-t_min+1,1e20)

        elif 6*t_min<=i<=6*t_max:
            for tau in range(t_min,t_max+1):
                if 5*t_min<=i-tau<=5*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i]
            temp=np.full(t_max-t_min+1,1e20)

        elif 7*t_min<=i<=7*t_max:
            for tau in range(t_min,t_max+1):
                if 6*t_min<=i-tau<=6*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i]
            temp=np.full(t_max-t_min+1,1e20)

        elif 8*t_min<=i<=8*t_max:
            for tau in range(t_min,t_max+1):
                if 7*t_min<=i-tau<=7*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i]
            temp=np.full(t_max-t_min+1,1e20)

        elif 9*t_min<=i<=9*t_max:
            for tau in range(t_min,t_max+1):
                if 8*t_min<=i-tau<=8*t_max:
                    #print temp,i,i-tau
                    temp[tau-t_min]= d[i-tau] + var_seg(inp[i-tau:i-1],u)
            d[i]=np.amin(temp)
            s[i]=np.append(s[i-(np.argmin(temp)+t_min)],int(i-(np.argmin(temp)+t_min)))
            #print i,s[i]
            temp=np.full(t_max-t_min+1,1e20)

        #print 'i', i, d[i]-d[i-1]
        if i == len(inp)  :
            factr=3

            temp2 =np.full(factr*t_max-factr*t_min+1,1e20)
            for last in range(factr*t_min,factr*t_max):
                temp2[last-factr*t_min] = d[last] + var_seg(inp[last:i-1],u)
                #print last, temp2[last-factr*t_min] ,d[last]
            
            #print np.argsort(temp2)+factr*t_min, np.sort(temp2)
            #print 'last',i,np.argmin(temp2)+4*t_min
            s[len(inp)]=np.append(s[np.argmin(temp2)+factr*t_min],np.argmin(temp2)+factr*t_min)
            #print s[i] 
            return s[i]
            break
        #if i ==len(inp):
        #    s[i]=np.append(s[i],i)
        #   return s[i]


"----------------function for min cost--------------------------------------"
def cost_fn(inp,u,s):
    sums=0
    for i in range(1,len(s)):
        sums=sums+var_seg(inp[s[i-1]:s[i]-1],u)
    return sums

"----------------update template -------------------------------------------"
def template(inp,seg):
    prev_cost = 1e20
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

"---------------iterations -------------------------------------------------"
def main(inp):
#    inp =der(inp)
    #u = np.zeros(t_min)
    u=inp[0:t_max-1]
    s = []
    cost=np.zeros(8)
    prev_cost =1e20
    for i in range(0,2):
        #print i
        s =segment(inp,u)
        #print cost_fn(inp,u,s)
        print s 
        u =template(inp,s)
        #print u, len(u)
        cost[i+1]= cost_fn(inp,u,s)
        print cost[i+1], len(u)
        #plt.plot(u)
        #plt.show()
        #plt.plot(inp)
        #for t in range(0,len(s)):
        #    plt.axvline(s[t],color='r',linestyle='--')
        #plt.show()
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
    #print freq,'Hz' 
    #s=np.delete(s,len(s)-1,0)
    #plt.plot(inp)
    #plt.plot(inp[s],s,'bo')
    #plt.show()


    
main(inp)
#s=[200,400,600,800]
#t=np.arange(0,1000,1)
#plt.plot(inp)
#plt.plot(s,inp[s],'gD')
#for i in range(0,len(s)):
#	plt.axvline(s[i],color='r',linestyle='--')
#plt.plot(inp[s],s,'bo')
#plt.show()
