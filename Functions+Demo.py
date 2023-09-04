import numpy as np
import pandas as pd
from scipy import signal
from scipy import stats
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import networkx.algorithms.community as nx_comm
%matplotlib qt

###################################################
#Define Functions
################################################### 
def normalize(vec):
    out=(vec-np.min(vec))/(np.max(vec)-np.min(vec))
    return(out)
###################################################
#Create Correlation-Matrices
################################################### 

def create_CorMatrix(df,plot=False):
    cor_df=pd.DataFrame(df) #assure that df is in pandas format
    cor_matrix=df.corr()**2 #cretae correlation matrix
    if plot==True: #define plot parameters
        ticks=np.linspace(0,len(cor_matrix)-1,5).astype("int")
        plt.imshow(cor_matrix,origin="lower",vmin=np.min(np.min(cor_matrix)),vmax=1)
        plt.title("Correlation Matrix",fontsize=14, fontweight="bold")
        plt.xlabel("Frequency")
        plt.ylabel("Frequency")
        plt.yticks(ticks,labels=cor_matrix.index[ticks])
        plt.xticks(ticks,labels=cor_matrix.index[ticks],rotation=45)
        cbar=plt.colorbar()
        cbar.ax.set_ylabel('R$^2$',fontsize=12, fontweight="bold",rotation=90)
        plt.tight_layout()
        plt.show()
    return cor_matrix

####################################################
#Find threshold for Correlation Matrix and create Adjacency Matrix
####################################################

def find_threshold(CorMatrix,variable_thresh=True,thresh=0,plot=False):
    if variable_thresh==True:
        cor_list=CorMatrix.to_numpy().flatten()
        cor_list=cor_list[cor_list!=1]
        g=np.append(cor_list,0)
        n = np.linspace(np.min(cor_list), np.max(cor_list), 100)
        nparam_density = stats.kde.gaussian_kde(g,bw_method="scott")
        dens_distribution = nparam_density(n)
        peak,_=signal.find_peaks(dens_distribution, prominence=0.1)
        low,_=signal.find_peaks(-dens_distribution)
        if len(low) > 0 and len(peak)>1:
            thresh_low=[low[-1]]
            threshold=n[low[-1]]
        else:
            e=dens_distribution[peak[-1]]/np.e
            x=np.where(dens_distribution==dens_distribution[np.abs(dens_distribution-e).argmin()])
            thresh_low=int(x[0])
            threshold=n[int(x[0])]
    if variable_thresh==False:
        threshold=thresh
    if plot==True:
        plt.plot(dens_distribution)
        plt.title("Distribution of Pearson Correlation Values",fontsize=14, fontweight="bold")
        plt.vlines(thresh_low,ymin=0,ymax=np.max(dens_distribution),color="red", label="Threshold = "+str(np.round(threshold,2)))
        ticks=np.linspace(0,99,10).astype("int")
        plt.xticks(ticks,labels=np.round(n[ticks],2))
        plt.xlabel("Pearson Correlation")
        plt.ylabel("Frequency")
        plt.legend(loc="upper left")
    return(threshold)


def create_AdjMatrix(CorMatrix,threshold,plot=False):
    cor_matrix=np.array(CorMatrix)
    np.fill_diagonal(cor_matrix, 0)
    cor_matrix=np.where(cor_matrix>=threshold,cor_matrix,0)
    cor_matrix=np.where(cor_matrix<threshold,cor_matrix,1)
    adj_matrix=cor_matrix
    adj_matrix=pd.DataFrame(cor_matrix,index=CorMatrix.index,columns=CorMatrix.index)
    if plot==True:
        ticks=np.linspace(0,len(cor_matrix)-1,5).astype("int")
        plt.imshow(adj_matrix,origin="lower",cmap="binary")
        plt.title("Adjacency Matrix",fontsize=14, fontweight="bold")
        plt.xlabel("Frequency")
        plt.ylabel("Frequency")
        plt.yticks(ticks,labels=CorMatrix.index[ticks])
        plt.xticks(ticks,labels=CorMatrix.index[ticks],rotation=45)
        plt.tight_layout()
        plt.show()
    return(adj_matrix)

###################################################
#Calculate Network Measures
###################################################    
def create_AssortDict(AdjMatrix):
    ind=AdjMatrix.index.values
    AssortDict={}       
    i=0
    for char in ind:
        if char.endswith("Hz")==True:
            n=char.strip("Hz")
        if char.endswith(" Hz")==True:
            n=char.strip(" Hz")
        else:
            n=char
        d={i:int(n)}
        AssortDict.update(d)
        i=i+1
    return(AssortDict)


def create_NetMeasures(AdjMatrix,AssortDict,plot_Com=False):
    adjmatr=np.array(AdjMatrix)
    net=nx.convert_matrix.from_numpy_matrix(adjmatr) #create network from AdjMatrix
    NetMeasures=pd.DataFrame(index=["Network Measures"],columns=["Assortativity","Average Shortest Path Length","Link Density","Modularity"])
    #assortativity
    nx.set_node_attributes(net, AssortDict, "Hz")
    try:
        NetMeasures["Assortativity"]=nx.numeric_assortativity_coefficient(net, "Hz")
    except:
        NetMeasures["Assortativity"]=np.nan
    #average shortest path length
    try:
        NetMeasures["Average Shortest Path Length"]=nx.average_shortest_path_length(net)
    except:
        NetMeasures["Average Shortest Path Length"]=np.nan
    #degree centrality
    Degree_Centrality=nx.degree_centrality(net)
    #link density
    NetMeasures["Link Density"]=nx.density(net)
    #modularity
    greedy_com=list(greedy_modularity_communities(net))
    NetMeasures["Modularity"]=nx_comm.modularity(net,greedy_com)
    if plot_Com==True:
        com_matrix=pd.DataFrame(0,columns=np.arange(len(AdjMatrix)),index=np.arange(len(AdjMatrix))) 
        c = list(greedy_modularity_communities(net))
        c_len=pd.DataFrame(index=np.arange(len(c)),columns=['size'])
        for l in range(0,len(c)):
            f=pd.DataFrame(l+1,columns=c[l],index=c[l]) #create dataframe for each cluster with respective freq bin, l=value set for community
            com_matrix=com_matrix.combine(f,np.add,fill_value=0) # combine           
            c_len.iloc[l,0]=(len(c[l]))
        ticks=np.linspace(0,len(AdjMatrix)-1,5).astype("int")
        plt.imshow(com_matrix,origin="lower",cmap="rainbow")
        plt.title("Greedy Communities (n="+str(len(c))+")",fontsize=14, fontweight="bold")
        plt.xlabel("Frequency")
        plt.ylabel("Frequency")
        plt.yticks(ticks,labels=AdjMatrix.index[ticks])
        plt.xticks(ticks,labels=AdjMatrix.index[ticks],rotation=45)
        plt.tight_layout()
        plt.show()           
    return(NetMeasures,Degree_Centrality)

###################################################################################################
###################################################################################################
                                        #Demo
###################################################################################################
###################################################################################################
                            
                            
###################################################
#Load Dataset containing power spectra information. 
###################################################
"""
For ideal output, the column names should represent the freqeuncy band and the index the audio recording ID.
Example:
---------------------------------------------------------------------
     ID      |     0 Hz     |  21 Hz   | ... |  22028 Hz | 22050 Hz |
---------------------------------------------------------------------
Recording_0 | 3.983402e-14 | 0.231118 | ... |  0.000421  | 0.000424 |
                             ........                            
Recording_n | 3.967565e-14 | 0.456618 | ... |  0.000601  | 0.000564 |
---------------------------------------------------------------------
"""
###################################################
#define directory where test_data.csv is located
path=""#enter path to test_data.csv
#load test_data.csv
df=pd.read_csv(path+"test_data.csv",index_col=0)

###################################################
#Create Correlation Matrix from Dataset of Power Spectra.
###################################################
CorMatrix=create_CorMatrix(df,plot=True)
"""
---------------------------------------------------------------------
create_CorMatrix(df,plot=False)
---------------------------------------------------------------------
Parameters:
df:               Dataframe of Power Spectra (for correct format see above).
plot:             If True, a plot of the Correlation Matrix is shown.
---------------------------------------------------------------------
Returns:
CorMatrix:        Correlation Matrix from Dataset of Power Spectra.
---------------------------------------------------------------------
"""
###################################################
#Find Threshold to define Adjacency Matrix.
###################################################
threshold=find_threshold(CorMatrix,plot=True)
"""
---------------------------------------------------------------------
find_threshold(CorMatrix,variable_thresh=True,thresh=0,plot=False)
---------------------------------------------------------------------
Parameters:
CorMatrix:         Correlation Matrix to threshold.
variable_thresh:   If True, the threhold will be found as defined in the paper.
thresh:            If variable_thresh=False, this value gives the threshhold instead.Must be between -1 and 1.
plot:              If True, a plot of the Pearson Correlation Value Distribution with the defined threshold is shown.
---------------------------------------------------------------------
Returns:
Threshold:         Float value of threshold.
---------------------------------------------------------------------
"""
###################################################
#Create Adjacency Matrix.
###################################################
AdjMatrix=create_AdjMatrix(CorMatrix,threshold,plot=True)
"""
---------------------------------------------------------------------
create_AdjMatrix(CorMatrix,threshold,plot=False)
---------------------------------------------------------------------
Parameters:
CorMatrix:         Correlation Matrix to threshold.
threshold:         Threshold that will be used to create Adjacency Matrix from the Correlation Matrix. All values below this value become zero, all other values become one.
plot:              If True, a plot of the Adjacency Matrix is shown.
---------------------------------------------------------------------
Returns:
AdjMatrix:        Pandas Dataframe of the thresholded Adjacency Matrix.
---------------------------------------------------------------------
"""
###################################################
#Create Assortativity Dictionary.
###################################################
AssortDict=create_AssortDict(AdjMatrix)
"""
---------------------------------------------------------------------
create_AssortDict(AdjMatrix)
---------------------------------------------------------------------
Parameters:
CorMatrix:         Adjacency Matrix or Correlation Matrix with the correctly formatted Hz values as index. This works for three Formats: XX Hz, XXHz, XX
---------------------------------------------------------------------
Returns:
AssortDict:        Dictionary of Hz Values for each node of AdjMatrix.
---------------------------------------------------------------------
"""
###################################################
#Calculate Network Measures.
###################################################
NetMeasures,Degree_Centrality=create_NetMeasures(AdjMatrix,AssortDict,plot_Com=True)
"""
---------------------------------------------------------------------
create_NetMeasures(AdjMatrix,AssortDict,plot_Com=False)
---------------------------------------------------------------------
Parameters:
AdjMatrix:         Adjacency Matrix to calculate Network Measures from.
AssortDict:        Dictionary of Hz Values to calculate Assortativity.
plot_Com:          If True, a plot of the Greedy Communities is shown. It should be noted, that this visualization is only feasible for Networks with high Assortativity.
---------------------------------------------------------------------
Returns:                
NetMeasures:       Pandas Dataframe with Assortativity, Average Shortest Path Length, Link Density, Modularity.
Degree_Centrality: Dictionary of Degree Centrality for each freqeuncy bin (i.e. node).
---------------------------------------------------------------------
"""
