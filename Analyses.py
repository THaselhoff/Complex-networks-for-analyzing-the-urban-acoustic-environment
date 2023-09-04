from scipy import signal						# version 1.5.2
import numpy as np 							# version 1.20.3
import pandas as pd 							# version 1.1.3
from scipy import stats							# version 1.5.2
import tqdm								# version 4.50.2
import networkx as nx							# version 2.5
from networkx.algorithms.community import greedy_modularity_communities	# version 2.5
import networkx.algorithms.community as nx_comm				# version 2.5
#%matplotlib qt


def normalize(vec):
    out=(vec-np.min(vec))/(np.max(vec)-np.min(vec))
    return(out)

###################################################
#Define AAD IDs for all recording devices
###################################################
AAD_ID=[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
###################################################
#Load Dataset of frequency spectra as described in 2.1 from each device and create Frequency Correlation Matrices (FCM)
################################################### 

for a in tqdm.tqdm(AAD_ID):
    r_list=[] #create list to store correlation matrices
    time=np.load("time.npy", allow_pickle=True) #read time information for each recording by device
    time=pd.DataFrame(time)
    time["Hour"]=time.iloc[:,0].dt.hour #generate variable with hour of day
    freq=np.load("signal.npy") #load frequency information for each recording by device (see 3.1)
    df=pd.DataFrame(freq[:-2,:]) #remove unneccessary meta information from f
    df_tf=pd.concat([time,df.T],axis=1) #concatenate time and frequency informaiton
    for i in tqdm.tqdm(range(0,24)): #create 24 Frequency Correlation Matrices by device, one for each hour of day
        df_hour=df_tf[df_tf.Hour==i]
        cor=df_hour.iloc[:,2:].corr()
        r_list.append(cor)
    np.save("FCM_"+str(a)+".npy",r_list) #save one file with 24 FCMs, seperately for each recording device

####################################################
#Define the threshold (see 3.2.2)
####################################################
peakandlow=pd.DataFrame(columns=["Device","Hour","Peak","Low"]) #create dataset to store peaks and lows, from which the thresholds are derived

kernel_dens_list=[] #create list to store kernel density estimations

#defining lists to store the respective information in the following loop
hour=[]
device=[]
lows=[]
peaks=[]
thresh_l=[]
thresh_p=[] 
z=0
for a in tqdm.tqdm(AAD_ID):
    time=np.load("time.npy",allow_pickle=True) #read time information for each recording by device
    z=0  
    kernel_dens=[]
    #Calculate peaks and lows for all FCMs by hour of day and per device
    for i in tqdm.tqdm(time):
        cor=i.flatten()**2
        g=cor[cor!=1] 
	g=np.append(cor,0)
        x = np.linspace(np.min(g), np.max(g), 100) #define sample from kernel density
        nparam_density = stats.kde.gaussian_kde(g,bw_method="scott") #estimate kernel density
        nparam_density = nparam_density(x) #sample from kernel density
        kernel_dens.append(nparam_density)
        nparam_density=normalize(nparam_density)
        f,_=signal.find_peaks(nparam_density, prominence=0.04) 
        l,_=signal.find_peaks(-nparam_density)
        thresh_p.append(f[-1])
        #define threshold using the criteria from 3.2.2 & Figure 2
        if len(l) > 0 and len(f)>1:
            thresh_l.append(l[-1])
        else:
            e=nparam_density[f[-1]]/np.e
            x=np.where(nparam_density==nparam_density[np.abs(nparam_density-e).argmin()])
            thresh_l.append(int(x[0]))
        hour.append(z)
        device.append(a)
        z=z+1    
    kernel_dens_list.append(kernel_dens)
  
peakandlow["Device"]=device    
peakandlow["Hour"]=hour    
peakandlow["Peak"]=thresh_p    
peakandlow["Low"]=thresh_l   

peakandlow.to_pickle("peak_and_low.pkl")
np.save("kernel_desnsities.npy",kernel_dens_list)   

####################################################
#Define Networks and create Adjacency Matrix (see 3.2.1, 3.2.2 and Figure 2)
####################################################
peakandlow=pd.read_pickle("peak_and_low.pkl")

for a in tqdm.tqdm(AAD_ID):
    adj=[] #create list to sotre Adjacency Matrices into
    z=0
    fcm=np.load("FCM_"+str(a)+".npy") #read 24 FCMs for each recording device
    t=peakandlow[peakandlow.Device==a]  #subsetting peakandlow, as only "Low" stores the threshold parameter
    for i in fcm:
        r=np.median(t.Low/100) #R-square thershold
        g=i**2 #square the FCM
        g=np.where(g>=r,g,0) #everything below r becomes 0
        g=np.where(g<r,g,1) #everything above r becomes 1
        np.fill_diagonal(g, 0) #remove the 1s caused by correlation of identical frequencies 
        adj.append(g)
        z=z+1
    np.save("AdjacencyMatrix_"+str(a)+".npy",adj) #save one file with 24 Adjacency Matrices, seperately for each recording device


###################################################
#Calculate Average Shortest Path Length (ASPL) (see 3.2.3)
###################################################

path_len_list=[]
for a in tqdm.tqdm(AAD_ID):
    adjmatr=np.load("AdjacencyMatrix_"+str(a)+".npy") #read 24 Adjacency Matrices for each recording device 
    path_len=[]
    for i in range(0,24):    
        net=nx.convert_matrix.from_numpy_matrix(adjmatr[i]) #define the graph from the Adjacency Matrix
        try:
            path_len.append(nx.average_shortest_path_length(net)) #calculate ASPL
        except:
            path_len.append(np.nan)
    path_len_list.append(path_len)
    
np.save("ASPL.npy",path_len_list) #save one file with all values of ASPL
















