from Bio import SeqIO
import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from sklearn.metrics import pairwise_distances
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""_This software takes aligned ribosomal gene sequences as input, sequencing error are masked by a cutomized
error rate and partially corrected during merging. The pairwise distance matrix between the reading are 
calculated. The indentical sequences are merged followed by merging between highly similar sequences. The outputs
including the pairwise distance after 0 distance merging , the merging outcome(ids and copynumbers) and aligned sequences in fasta file._
Args:
     input_file (filepath): The fasta file
Output_files:
    1. ceratain readings: pandas csv of ceratain readings:with index,seqID and copy number(C_Num)
    2. ceratain+uncertain readings: pandas csv of ceratain+uncertain readings, with index,seqID and copy number(C_Num)
    3. final_merged_certain: fasta: sequences of ceratain readings
    4. final_merged_certananduncertain: fasta: sequences of ceratain and uncertain readings
    5. Counter_Pd: pandas Counter csv: all readings with SeqID and C_Num after 0 distance merging
    6. fasta: sequences after 0 distance merging
    7. Mdis: pairwise distance matrix

"""
def CreateMatrix(input_file):
    """This function creates a matrix from a given fasta file. 

    Args:
        input_file (filepath): The fasta file
    
    return: M(2d array of stringsy): matrix with M_rows of rows and M_cols of columns.
    """
   
    records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    ids = list(SeqIO.to_dict(SeqIO.parse(input_file, "fasta")).keys())
    array =[
        y
        for i in range(0,len(ids))
            for y in records[ids[i]].seq]
        
    m = np.array(array)
    M = np.reshape(m, (-1, len(records[ids[1]].seq)))
    M_rows, M_cols = M.shape 
    
    return M,M_rows,M_cols
    
def MarkTerminalGaps(M,M_rows, M_cols):
    """_The function Mark the terminal gaps, gaps'-' within the terminal gaps are masked by '+'.
        Defining terminal gaps: 
        1:from the begining or the end of the readings;
        2: at least 3 in a row;
        3: terminated with the first non-gap character;

    Args:
        M (2d array of stringsy):Sequence matrix, product of CreateMatrix(input_file)
        M_rows (int): :rows of array of strings, product of CreateMatrix(input_file)
        M_cols (int):columns of array of strings, product of CreateMatrix(input_file)
    Returns:
        M (array of strings):Sequence matrix
        M_rows (int ): :rows of array of strings
        M_cols (int):columns of array of strings
        MissInfoRowCoords(int array): coordinates of all the readings with terminal gaps
    """
    ##find all gaps
    gapRowCoords,gapColCoords = np.where(M == '-')
    gapCoords=np.vstack((gapRowCoords,gapColCoords)).T


    #fw:the rows with 1st 3 chars are gaps
    MissInfoRowCoordsF=[
        (i)
    for i,j in gapCoords
        if j==0 and M[i,(j+1)]=='-' and M[i,(j+2)]=='-'
]    

    #Bw:the rows with last 3 chars are gaps
    MissInfoRowCoordsB=[
        (i)
        for i,j in gapCoords
            if j==M_cols-1 and M[i,(j-2)]=='-'and M[i,(j-3)]=='-'
] 

    #conbine Fw and Bw
    MissInfoRowCoords=list(set(MissInfoRowCoordsF + MissInfoRowCoordsB))
        # mask rows with terminal gaps 
    for i in MissInfoRowCoordsF:
        for j in range(0,M_cols):
            if M[i,j]== '-':
                M[i,j] = '+'
            if M[i,j+1]!= '-':
                break
    for i in MissInfoRowCoordsB:
        for j in reversed(range(0,M_cols)):
            if M[i,j]== '-':
                M[i,j] = '+'
            if M[i,j-1]!= '-':
                break
   
    M_rows,M_cols = M.shape
    print(len(MissInfoRowCoords)," rows reads with terminal gaps have been Masked.")
    return M,M_rows,M_cols,MissInfoRowCoords

def ErrorMasking(M,M_rows,M_cols,T):
    """This function masks the error reads by analyzing the majority of readings on the same postions, error rate is  
    one of the Input value.

    Args:
        M (2d array of stringsy):Sequence matrix
        M_rows (int ): :rows of array of strings
        M_cols (int):columns of array of strings
        T (Int): customer defined error rate
    return: 
        n/a: the error masking is done onto the matrix of readings.
    """
    for x in range(0,M_cols):
        ele_array = np.unique(M[:,x])
    
        if len(ele_array) == 1:
                pass
        else:
               for i in ele_array:
                ele_per = ((M[:,x]==i).sum(axis=0))/M_rows
           
                if ele_per > T:
                    pass
                elif ele_per <= T and ((M[:,x]=='-').sum(axis=0))>(M_rows*(9/10)):
                    for j in (np.where(M[:,x]==i)):
                        M[j,x]='+'
                else:
                    for j in (np.where(M[:,x]==i)):
                        M[j,x]='n'
            

def StringMatrixToNumber (M1): 
    
    """_Transfer string array in to Int array. Each base and gaps are represented by an int_:
        1: a
        2: c
        3: t
        4: g
        n: 0
        -: -2
        +: -1

    Args:
        M1 (2d array of strings):

    Returns:
        M2(2d array of int):
    """
    
    M2=M1
    M2[M2 == 'a']="1"
    M2[M2 == 'c']="2"
    M2[M2 == 't']="3"
    M2[M2 == 'g']="4"
    M2[M2 == 'n']="0"
    M2[M2 == '-']="5"#-2
    M2[M2 == '+']="6"#-1
    Int = np.vectorize(lambda x: int(x))
    M2 = Int(M1)   
    M2[M2 == 5]=-2
    M2[M2 == 6]=-1
    print("First three position of number Matrix:",M2[0,0],M2[0,1],M2[0,2],"Datatype:",M2.dtype)
    return M2
    
def CalculatePairWiseDistance_sk(u,v):
    """_customizing calculation formula for sklearn pairwise distance matrixfunction_

    Args:
        u (1d array of strings): 1st row of seq reading
        v (1d array of strings): the row of seq reading to pair with 'u' for distance calculation

    Returns:
        a+b(int):the distanct between u and v
    """    
    a = np.bitwise_and(np.bitwise_and(u > 0,v > 0),(u!=v)).sum(axis=0)
    b = np.bitwise_and(u == -2,v>0).sum(axis=0) + np.bitwise_and(v == -2,u>0).sum(axis=0)
    return a+b

def CalculatePairWiseDistance(M2):
    """_This function generate a counter that records the copy number of the number matrix 
        after indentical merging. When two sequences are indentical, the copyNumber 
        of the first template reading increased by 1 and the copynumber of the second seq reading is marked as 0_
        RepeatSystem: In order to avoid repeative calculation onto the indentical sequences. A repeat switch 
        has been implemented, 
        'n' and '+' correction: when indentical merging happens. The 'n's and '+' on the template seq reads are tried to be 
        corrected by the base from the same position of its merging template. If a 'non-n' base appears on the same position
        'n's and '+'s are corrected by non-n bases
    Args:
        M2 (2d array of int): 

    Returns:
        counter(int array):a counter that record the copy number of the number matrix
        M2(int matrix): int matrix, 'n' and '+' corrected
        # Mdis(fullscale): pairwise distance matrix of all readings including 0 distance, can be returned if needed
        # MergeDicts_Pd#Mdis(dictionary of int): a dictionary of the merging path can be calculated and exported in needed
    """    
    
    M1=M2
    M_rows,M_cols = M2.shape
    # distance matrix with M_rows* M_rows shape
    Mdis = np.full((M_rows,M_rows),-1)
    repeat = np.full(M_rows,1,dtype=bool)
    counter = np.full(M_rows,1,dtype=int)
    #MergeDicts_Pd = pd.DataFrame( {'SeqID_i':[],  'SeqID_j':[]})
    nmcounter =0
    for i in range(0,M_rows-1):
        if repeat[i]==False:
                continue
        else:
            for j in range(i+1,M_rows):
                if repeat[j]==False:
                     continue
                else:
                    a = np.bitwise_and(np.bitwise_and(M1[i,]>1,M1[j,]>1),(M1[i,]!=M1[j,])).sum(axis=0)
                    b = np.bitwise_and(M1[i,]== -2,M1[j,]>0).sum(axis=0) + np.bitwise_and(M1[j,]== -2,M1[i,]>0).sum(axis=0)
                    Mdis[i,j]= a+b  
                    if a+b == 0:
                        counter[i]+=1
                        counter[j]-=1
                        repeat[j]= False
                        for x in range(0, len(M1[i,])):
                            if M1[i,x] == 0 or M1[i,x] == -1:
                                M1[i,x] == M1[j,x]
                                nmcounter += 1

    #1.1 If we want to see How Merge happens unhash the code below
    #df2 = pd.DataFrame([[i,j]], columns=['SeqID_i','SeqID_j'])
    # MergeDicts_Pd = pd.concat([df2, MergeDicts_Pd])
                      #  else:
    print(nmcounter,'of ns and gaps has been corrected')             
    return counter,M2 #MergeDicts_Pd#Mdis


def reduction_Matrix(M2,counter):
    """_This function remove the merged readings fromt he 2d array of the readings in int dtype based on the counter.  
    The readings with postive copynumbers were returned as output_

    Args:
        M2 (_2d array of int_): 
        counter (_1d array of int_): output of the function'CalculatePairWiseDistance(M2)'

    Returns:
        array2(2d array of int): M2 after removing merged reading.
    """    
    M_rows,M_cols = M2.shape
    dict1={}
    for i,j in enumerate(counter):
        dict1[i]=j
    array2 = np.zeros((len(dict1.keys()),M_cols))
    for i in range(0,len(dict1.keys())):
        if dict1[i]==0:
            continue
        else:
            array2[i]=M2[list(dict1.keys())[i]]
    return array2

def Export(Mdis):
    """function that exporting pairwise distance matrix into csv file, more csvs such as merging path can be exported if necssary.

    Args:
        Mdis (2d array of int), 
        
    """    
    from numpy import savetxt
    savetxt('Dist_Matrix_v2.0_reduction_1e2.csv',Mdis,fmt='%i',delimiter=',')
    #savetxt('Counter1e2.csv',counter,fmt='%i',delimiter=',')
    #savetxt('Merge_list.csv',Merge_list,fmt='%i',delimiter=',')#   continue


def NumberMatrixToString (reducedarray):
    """_This function transfer 2d array of int back to 2d array of strings_

    Args:
        reducedarray (2d array of int): 

    Returns:
        array_rev(2d array of srtings): 
    """    
    rows,cols = reducedarray.shape
    array_rev = np.full((rows,cols),'a')
    for i in range(0,rows):
        for j in range(0,cols):
            if reducedarray[i,j]==1.0:
                array_rev[i,j]='a'
            elif reducedarray[i,j]==2.0:
                array_rev[i,j]='c'
            elif reducedarray[i,j]==3.0:
                array_rev[i,j]='t'
            elif reducedarray[i,j]==4.0:
                array_rev[i,j]='g'
            elif reducedarray[i,j]==0.0:
                array_rev[i,j]='n'
            else:
                array_rev[i,j]= '-'
    return array_rev
    
def StringMatrixToSeq(array_rev):
    """_This funciton removes the ',' of 2d array of string, reduce the array dimention to 1d _
        input: ['a','c','t','g'], output:['actg'] 
    Args:
        array_rev (2d array of strings)

    Returns:
        seq_array(1d array of strings): 
    """    
    rows,cols = array_rev.shape
    separator = ''
    seq_array = np.array([
            separator.join(array_rev[i,])
            for i in range(0,rows)
        ])
    return seq_array

def FastaRecords(seq_array):
    """_This function write 1d array of strings into SeqRecords of SeqIO_, and expor a fasta file

    Args:
        seq_array (1d array of strings)
    """    
    rows, = seq_array.shape
    separator = ''
    a = np.array([seq_array[i].item() for i in range(0,rows)])
    b = [Seq(i) for i in a]
    records = [
           SeqRecord(Seq(seq), id = str(index), description = "") 
           for index,seq in enumerate(b) 
    ]
    SeqIO.write(records, 'Merged_readings' ,"fasta")
    return 
# a counter df with seqID
def GenerateCounterDataFrame(counter):
    """ This function generates a counter dataframe for each reading and with seqID and copy number(C_Num).

    Args:
        counter (1d array of int): Counter of copyNumbers

    Returns:
        Counter_Pd(pandas): DataFrame with Columns: CopyNumber:C_Num and SeqIDs
    """   
  
    Counnter_Dict = {}
     ### NonZeroArray: Counter without 0 copyNumbers

    NonZeroDict = {}
    
    for i,j in enumerate(counter):
        Counnter_Dict[i] = j
    
    for i in range(0,len(Counnter_Dict.values())):
        if counter [i]!=0:
            NonZeroDict[i] = Counnter_Dict[i]
    
    NonZeroArray = np.array(list(NonZeroDict.values()))
    
    Counter_Pd = pd.DataFrame(NonZeroArray)
    Counter_Pd = Counter_Pd.rename(columns = { 0:'C_Num'})
    Counter_Pd['SeqID'] = range(0,len(Counter_Pd))
    Counter_Pd.to_csv('Counter_Pd.csv',sep=',')
    return Counter_Pd

def DeleteIsolateSingleReadings(Mdis,Counter_Pd,T):
    """Delete the isolated single readings that are ditanct with any other readings.
        If a reading's at least 5 step away from any of the other readings, then is will be removed

    Args:
        Mdis (2d array of pairwise_distances):
        Counter_Pd (dataframe of seq readings):
        T (int): distanct threshold that defines distant
    Returns:
        Counter_Pd_removeDist(dataframe of seq reading): Counter_Pd removed ditant readigns
    """    
    str = f"MdisM{T}" 
    exec("%s = Mdis<T " % (str))
    DistList, = np.where(eval(str).sum(axis=0)==1)
    DistList_CNumblowT = np.array([
        i
    for i in DistList
        if int(Counter_Pd.loc[Counter_Pd["SeqID"]==i,"C_Num"]) == 1
]) 
    Counter_Pd_removeDist = Counter_Pd.drop(DistList_CNumblowT)
    return Counter_Pd_removeDist
def OneCopyReadingList(Counter_Pd_Unique_RemoveUniqueMerge_removeDist):
    """_The merging logic after 0 dist merging is merging the 1 copy readings to the multiple copy one. This function locate all the 1 copy readings and 
    organize then into a pandas dataframe_

    Args:
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist (dataframe of seq reading): Counter_Pd removed ditant readigns

    Returns:
        OneCopyDf(pandas dataframe): dataframe of 1 copy readings with SeqID and C_Num
        OneCopyDfSeqID(numpy 1d array): array of the 1 copy seqIDs
    """    
    OneCopyDf = Counter_Pd_Unique_RemoveUniqueMerge_removeDist[Counter_Pd_Unique_RemoveUniqueMerge_removeDist["C_Num"]== 1]
    OneCopyDfSeqID = np.array(list(OneCopyDf["SeqID"]))
    return OneCopyDf,OneCopyDfSeqID
def CoordsOfOnedistInDistMatrix(Mdis):
    """This function generate a dataframe of seqIDs which are 1 step away from each other.

    Args:
        Mdis (2d array of pairwise_distances)
    Returns:
        OneDist_pd(dataframe of 2 cols seqIDs): dataframe is representing the connections in columns of seqIDi and SeqIDj.
    """    
    def unique_rows(a):
        """_A function that take an array tuples as input and remove the duplicated tuples within the tuple array_. 
        Reason of implementation:distance coords of the shows up repeatly from the squared pairwise distance matix. In order to avoid 
        repeatly calculation of the same pairs.

        Args:
            a (_1d array of tulple ): 

        Returns:
              1d array of unique tuples.
        """        
        a = np.ascontiguousarray(a)
        unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
        return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    
    r,c = np.where(Mdis == 1)
    XCoords=np.vstack((r,c)).T
    XCoords = unique_rows(np.sort(XCoords))
    OneDist_pd = pd.DataFrame(XCoords).rename(columns = { 0:'Seq_ID_i',1:'Seq_ID_j'})
    return OneDist_pd

def FindUniqueOneDistReading(OneDist_pd,Counter_Pd):
    """_This function find readings that uniquely link to each other and merge the readings by adding 1 to the template's copynumber, the ones which has merged
    to other readings are recorded into a droplist waiting for dropping._

    Args:
        OneDist_pd (pandas dataframe ): _dataframe of 1dist connection_
        Counter_Pd (pandas dataframe ): _datafrom with seqIDs and C_Num_

    Returns:
        droplist2_array(1d np array):seqIDs that have been merged into other readings
        Counter_Pd(pandas dataframe):_datafrom with seqIDs and C_Num, updated C_Num from merging_
    """    
    droplist2 = {}
    for i in range(0,len(Counter_Pd)):
        for j in OneDist_pd.loc[OneDist_pd["Seq_ID_i"]==i,"Seq_ID_j"]:
            if len(OneDist_pd.loc[OneDist_pd["Seq_ID_i"]==j,"Seq_ID_j"]) == 1:
                droplist2[j]= j
                Counter_Pd.loc[Counter_Pd['SeqID'] == i,['C_Num']] += int(Counter_Pd.loc[Counter_Pd['SeqID']==j,'C_Num'])
                droplist2_array = np.array(list(droplist2.keys()))
    return droplist2_array,Counter_Pd
def DefineUncertainReadings(Counter_Pd_Unique_RemoveUniqueMerge_removeDist,T1,T2):
    """_From Counter_Pd this function isolate the uncertain readings of interested:
        define uncertain readings of interested: the readings with C_Num between threshold T1 and T2. From the dataframe after Uniquemerging.

    Args:
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist (pandas datafrom with seqIDs and C_Num): datafrom with seqIDs and C_Num, updated C_Num from merging
        T1 (int): _lower threshold for uncertain readings of interested_
        T2 (in): _upper threshold for uncertain readings of interested_

    Returns:
        UncertainReadings(dataframe of Uncertain readings of Interest): Counter_Pd removed ditant readigns
        UncertainReadingsID(1d array of int):(array of seqIDs of uncertainreadings
    """    
    #Input: T1: lower bond of Uncertainreadings, T2: upperBond of UncertainReadings
    UncertainReadings = (Counter_Pd_Unique_RemoveUniqueMerge_removeDist[Counter_Pd_Unique_RemoveUniqueMerge_removeDist["C_Num"] >T1]).loc[(
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist[Counter_Pd_Unique_RemoveUniqueMerge_removeDist["C_Num"] >T1])['C_Num']<T2]
    UncertainReadingsID = np.array(list(UncertainReadings['SeqID']))
    return UncertainReadings,UncertainReadingsID
def Split_Merge_1copy_to_uncertain(UncertainReadings,UncertainReadingsID,OneCopyDfSeqID,Mdis,Counter_Pd,T):
    """_Split merging readings with 1 copies to uncertain readings_

    Args:
        UncertainReadings (_pdDataframe_): 
        UncertainReadingsID (_np1darray_): selected id array
        OneCopyDfSeqID (_np1darray_): selected id array
        Mdis (_2d array of int_): distance matrix
        Counter_Pd (_pdDataframe_): dataframe with SeqID and C_Num
        T (int): distance of of selection

    Returns:
       UncertainReadings(dataframe)
       droplist_PartiallyMergedToUncertain_array(1d array of int)
    """    
    Mdis_Pd = pd.DataFrame(Mdis)
    droplist_PartiallyMergedToUncertain = {}
    
    for i in UncertainReadingsID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                droplist_PartiallyMergedToUncertain[j] = j
                UncertainReadings.loc[UncertainReadings['SeqID'] == i,'C_Num'] +=  float(UncertainReadings.loc[UncertainReadings['SeqID'] == i,'C_Num']) /sum([
                float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                ])
                droplist_PartiallyMergedToUncertain_array = np.array(list(droplist_PartiallyMergedToUncertain.keys()))
    try:
        return UncertainReadings,droplist_PartiallyMergedToUncertain_array
    except:
        UncertainReadings,droplist_PartiallyMergedToUncertain_array = [],[]
        print('no seq can be merged to uncertains')
        return UncertainReadings,droplist_PartiallyMergedToUncertain_array
def DefineCertainlist(Counter_Pd_Unique_RemoveUniqueMerge_removeDist,T1):
    """_Define Certain readings list of interest. 
    Certain Reading: readings with C_Num that above T1_

    Args:
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist (pandas datafrom with seqIDs and C_Num): datafrom with seqIDs and C_Num, updated C_Num from merging
        T1 (int): _upper threshold for uncertain readings of interested_

    Returns:
        Certainlist(dataframe)
        CertainlistID(1d array of int)
    """    
  
    Certainlist = Counter_Pd_Unique_RemoveUniqueMerge_removeDist[Counter_Pd_Unique_RemoveUniqueMerge_removeDist["C_Num"] >=T1] 
    CertainlistID = np.array(list(Certainlist['SeqID']))
    return Certainlist,CertainlistID
def Split_Merge_1copy_to_certain(Certainlist,CertainlistID,OneCopyDfSeqID,Mdis,Counter_Pd,T):
    """Similar to Split_Merge_1copy_to_uncertain
    """    
    PartialMergedToCertain = {}
    Mdis_Pd = pd.DataFrame(Mdis) 
    for i in CertainlistID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                PartialMergedToCertain[j] = j
                Certainlist.loc[Certainlist['SeqID'] == i,'C_Num'] += float(Certainlist.loc[Certainlist['SeqID'] == i,'C_Num']) /sum([
                float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                ])
                droplist_PartialMergedToCertain_array = np.array(list(PartialMergedToCertain.keys()))
    try:
        return Certainlist,droplist_PartialMergedToCertain_array
    except:
        Certainlist,droplist_PartialMergedToCertain_array = [],[]
        print('no seq can be merged to certains')
        return Certainlist,droplist_PartialMergedToCertain_array
def OneCopyreadingNotJoinedMerging(OneCopyDf,OneCopyDfSeqID,droplist_PartiallyMergedToUncertain_array, droplist_PartialMergedToCertain_array):
    """_Output a dataframe  and an array list of 1 copy readings that haven't join any 1 dist merging_

    Args:
        OneCopyDf (dataframe): 
        OneCopyDfSeqID (1d array of seqIDs): 
        droplist_PartiallyMergedToUncertain_array (1d array of seqIDs): 
        droplist_PartialMergedToCertain_array (1d array of seqIDs): 

    Returns:
       OneCopyDfNotInvolved(dataframe)
       OneCopyReadingsNotInvolved(1d array of seqIDs)
    """    
    a = np.concatenate((droplist_PartiallyMergedToUncertain_array, droplist_PartialMergedToCertain_array), axis=None)
    ReadingsInvolved = np.unique(a)
    OneCopyReadingsNotInvolved = OneCopyDfSeqID[~np.isin(OneCopyDfSeqID, ReadingsInvolved)]
    OneCopyDfNotInvolved = OneCopyDf.drop(ReadingsInvolved)
    return OneCopyDfNotInvolved,OneCopyReadingsNotInvolved
def DataFrameFromOneDistMerging(Certainlist,UncertainReadings,OneCopyDfNotInvolved):
    """_DataFrame after 1 dist merging_

    Args:
        Certainlist (_dataframe_): 
        UncertainReadings (_dataframe_): 
        OneCopyDfNotInvolved (_dataframe_): 

    Returns:
        OneDistMergingDf(dataframe): DataFrame with Columns: CopyNumber:C_Num and SeqIDs
    """    
    OneCopyDfNotInvolved["C_Num"] = pd.to_numeric(OneCopyDfNotInvolved["C_Num"])
    a = Certainlist.merge(UncertainReadings,how='outer')
    OneCopyDfNotInvolved['C_Num'] = OneCopyDfNotInvolved['C_Num'].astype('float64')
    OneDistMergingDf = a.merge(OneCopyDfNotInvolved,how='outer')
    return OneDistMergingDf
    
    """_The following function are optimized merging functions for 2 and 3 distance merging, following the same logic as 1 dist merging. 
    starting and ending points marked
    """
### function for 2-3 dist merging starts here
def CoordsOfTdistInDistMatrix(Mdis,T,arrayofInterest):
    
    def unique_rows(a):
        a = np.ascontiguousarray(a)
        unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
        return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    
    r,c = np.where(Mdis == T)
    r1 = np.array([
        i
        for i in r
            if i in arrayofInterest 
    ])
    c1 = np.array([
        i
        for i in c
            if i in arrayofInterest 
    ])
    XCoords=np.vstack((r1,c1)).T
    XCoords = unique_rows(np.sort(XCoords))
    TDist_pd = pd.DataFrame(XCoords).rename(columns = { 0:'Seq_ID_i',1:'Seq_ID_j'})
    return TDist_pd
### Did not inheret funtion due to datatype: One dist was int, Two need float. may be able to find a way to sync both. try update in V3.1


def FindUniqueTwoDistReading(TwoDist_pd,OneDistMergingDf):
    droplist_2distUnique = {}
    for i in range(0,len(OneDistMergingDf)):
        for j in TwoDist_pd.loc[TwoDist_pd["Seq_ID_i"]==i,"Seq_ID_j"]:
            if len(TwoDist_pd.loc[TwoDist_pd["Seq_ID_i"]==j,"Seq_ID_j"]) == 1:
                droplist_2distUnique[j]= j
                OneDistMergingDf.loc[OneDistMergingDf['SeqID'] == i,['C_Num']] += float(OneDistMergingDf.loc[OneDistMergingDf['SeqID']==j,'C_Num'])
                droplist_2distUnique_array = np.array(list(droplist_2distUnique.keys()))
    return droplist_2distUnique_array,OneDistMergingDf
def d2Split_Merge_1copy_to_uncertain(UncertainReadings,UncertainReadingsID,OneCopyDfSeqID,Mdis,arrayofInterest,Counter_Pd,T):
    Mdis_Pd = pd.DataFrame(Mdis)
    droplist_PartiallyMergedToUncertain = {}
    
    for i in UncertainReadingsID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                droplist_PartiallyMergedToUncertain[j] = j
                try:
                    UncertainReadings.loc[UncertainReadings['SeqID'] == i,'C_Num'] +=  float(UncertainReadings.loc[UncertainReadings['SeqID'] == i,'C_Num']) /sum([
                    float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                    for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                        if k in arrayofInterest
                    ])
                except:
                    UncertainReadings.loc[UncertainReadings['SeqID'] == i,'C_Num'] += 0
                
                droplist_PartiallyMergedToUncertain_array = np.array(list(droplist_PartiallyMergedToUncertain.keys()))
    try:
        return UncertainReadings,droplist_PartiallyMergedToUncertain_array
    except:
        UncertainReadings,droplist_PartiallyMergedToUncertain_array = [],[]
        print('no seq can be merged to uncertains')
        return UncertainReadings,droplist_PartiallyMergedToUncertain_array
def d2Split_Merge_1copy_to_certain(Certainlist,CertainlistID,OneCopyDfSeqID,Mdis,Counter_Pd,arrayofInterest,T):
    PartialMergedToCertain = {}
    Mdis_Pd = pd.DataFrame(Mdis) 
    for i in CertainlistID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                PartialMergedToCertain[j] = j
                try:
                    Certainlist.loc[Certainlist['SeqID'] == i,'C_Num'] += float(Certainlist.loc[Certainlist['SeqID'] == i,'C_Num']) /sum([
                    float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                        for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                            if k in arrayofInterest
                    ])
                except:
                    Certainlist.loc[Certainlist['SeqID'] == i,'C_Num'] += 0
                droplist_PartialMergedToCertain_array = np.array(list(PartialMergedToCertain.keys()))
    try:
        return Certainlist,droplist_PartialMergedToCertain_array
    except:
        Certainlist,droplist_PartialMergedToCertain_array = [],[]
        print('no seq can be merged to certains')
        return Certainlist,droplist_PartialMergedToCertain_array
### function for 2-3 dist merging ends here


def MaskNwithXDist_pd(XDist_pd,array_rev,CertainAndUncertainseqIDarray):
    """_This function is able to find the 'n and the terminalgap'_'s within 1 distance connect frame and correct 'n's in the merging template readings by
    the readings merging to it.

    Args:
        XDist_pd (dataframe): 
        array_rev (1d array of strings):
        CertainAndUncertainseqIDarray (1d array of seqIDs):

    Returns:
        array_rev (1d array of strings):n's are corrected
    """    
    r,c  = array_rev.shape
    Ncounter = 0
    for x in range(0,len(XDist_pd)):
        for y in range(0,c):
            if XDist_pd.iloc[x]['Seq_ID_i'] in CertainAndUncertainseqIDarray:
                if array_rev[XDist_pd.iloc[x]['Seq_ID_i'],y] == 'n' or array_rev[XDist_pd.iloc[x]['Seq_ID_i'],y] == '-':#update should be done in 3.1
                    array_rev[XDist_pd.iloc[x]['Seq_ID_i'],y] = array_rev[(XDist_pd.iloc[x]['Seq_ID_j']),y]
                    Ncounter +=1
    print(Ncounter,'n and gaps has been corrected')
    return array_rev
        
def FastaRecords_Final_CandUnC(MergedSeqRecords_array):
    """_This function is able to out a Fasta file from cerainanduncertain list of data of interested _

    Args:
        MergedSeqRecords_array (arrayofMegedseq): 
    """    
    rows, = MergedSeqRecords_array.shape
    separator = ''
    a = np.array([MergedSeqRecords_array[i].item() for i in range(0,rows)])
    b = [Seq(i) for i in a]
    records = [
           SeqRecord(Seq(seq), id = str(index), description = "") 
           for index,seq in enumerate(b) 
    ]
    SeqIO.write(records, 'Final_Merged_CertainandUncertain' ,"fasta")
    return 
def FastaRecords_Final_Certain(MergedSeqRecords_array):
    """_This function is able to out a Fasta file from cerain list of data of interested _

    Args:
        MergedSeqRecords_array (arrayofMegedseq): 
    """    
    rows, = MergedSeqRecords_array.shape
    separator = ''
    a = np.array([MergedSeqRecords_array[i].item() for i in range(0,rows)])
    b = [Seq(i) for i in a]
    records = [
           SeqRecord(Seq(seq), id = str(index), description = "") 
           for index,seq in enumerate(b) 
    ]
    SeqIO.write(records, 'Final_Merged_Certain' ,"fasta")
    return 
#def FindSeqIDandCNum(seqIDarray2,CNumarray2,Index):
    """_This function is able to quickly find the seqID and C_Num by index, muted in python file, can be used in JupiterNotebook_

    Args:
        MergedSeqRecords_array (arrayofMegedseq): 
    """      
 #   try:
  #      print('SeqID of','Index',f'{Index}',',is',seqIDarray2[Index],',copy number is',CNumarray2[Index])
   # except:
    #    print('Index is out of range')
def ExportMergeFinal(certaindf,ThreeDistMergingDf):
    """_Export final dataframes_

    Args:
        certaindf (dataframe): 
        ThreeDistMergingDf (dataframe): 
    """    
    certaindf.to_csv('CertainReadings.csv',sep=',')
    a = ThreeDistMergingDf.loc[ThreeDistMergingDf['C_Num']!=1,]
    a.to_csv('CertainandUncertainReadings.csv',sep=',')
     
def main():
       
    import timeit
    start = timeit.default_timer()
    input_file = sys.argv[1]  
    #read data
    M,M_rows,M_cols = CreateMatrix(input_file)
    #marktermmincal gap and mark the terminal gasp as '+'
    M,M_rows,M_cols,MissInfoRowCoords = MarkTerminalGaps(M,M_rows, M_cols)
    #sort terminal gaps indexes
    MissInfoRowCoords = np.sort(MissInfoRowCoords)
    # delete terminal gaps from the main matrix
    M_WithOutTerminalGapReadings = np.delete(M,(MissInfoRowCoords),axis=0)
    # add the terminal gaps to the end fo the main matrix
    ReMergeM = np.concatenate((M_WithOutTerminalGapReadings,M[MissInfoRowCoords]),axis=0)
    M_rows,M_cols=ReMergeM.shape
    # Mask sequencing erros
    ErrorMasking(ReMergeM,M_rows,M_cols,0.01)
    # StringMatrixToNumber
    M_Num = StringMatrixToNumber (ReMergeM)
    # Counter after merging the 0 distance readings and correct 'n's while 0 dist merging
    counter,M2 = CalculatePairWiseDistance(M_Num)
    # remove the merged readings from the main matrix
    array2 = reduction_Matrix(M_Num,counter)
    reducedarray = array2[~np.all(array2 == 0, axis=1)]
    # NumberMatrixToString 
    array_rev = NumberMatrixToString (reducedarray) 
    # pairwise_distances without 0 distance
    Mdis = pairwise_distances(reducedarray,metric = CalculatePairWiseDistance_sk)
    # export pairwise_distances metrix in csv format
    Export(Mdis)
    # string matrix to seq_array, reduce array dimention for fasta file
    seq_array = StringMatrixToSeq(array_rev)
    # export fasta file of seq readings after 0 distance merging
    FastaRecords(seq_array)
    # pandas df of seqID and CopyNumber
    Counter_Pd = GenerateCounterDataFrame(counter)
    # readings with 1 distance between them
    OneDist_pd=  CoordsOfOnedistInDistMatrix(Mdis)
    # Find the FindUniqueOneDistReading dataframe
    droplist2_array,Counter_Pd_Unique = FindUniqueOneDistReading(OneDist_pd,Counter_Pd)
    # Delete Isolate Single Readings
    Counter_Pd_removeDist = DeleteIsolateSingleReadings(Mdis,Counter_Pd_Unique,5)
    # Removed the reading appear in droplist after 1 dist unique merging
    Counter_Pd_Unique_RemoveUniqueMerge_removeDist = Counter_Pd_removeDist.drop(droplist2_array)
    #1 copy readings that are going to merge to other readings
    OneCopyDf,OneCopyDfSeqID = OneCopyReadingList(Counter_Pd_Unique_RemoveUniqueMerge_removeDist)
    #define uncerting readings
    UncertainReadings,UncertainReadingsID = DefineUncertainReadings(Counter_Pd_Unique_RemoveUniqueMerge_removeDist,3,10)
    #split merging 1copy readings to uncerting readings
    UncertainReadings,droplist_PartiallyMergedToUncertain_array = Split_Merge_1copy_to_uncertain(UncertainReadings,UncertainReadingsID,OneCopyDfSeqID,Mdis,Counter_Pd,1)
    #define certainReadings
    Certainlist,CertainlistID = DefineCertainlist(Counter_Pd_Unique_RemoveUniqueMerge_removeDist,10)
    #split merging 1 copy readings to certainReadings.
    Certainlist,droplist_PartialMergedToCertain_array= Split_Merge_1copy_to_certain(Certainlist,CertainlistID,OneCopyDfSeqID,Mdis,Counter_Pd,1)
    # find 1 copy readings that did not join the 1 distance merging of unique or to certain or uncertain.
    OneCopyDfNotInvolved,OneCopyReadingsNotInvolved =  OneCopyreadingNotJoinedMerging(OneCopyDf,OneCopyDfSeqID,droplist_PartiallyMergedToUncertain_array, droplist_PartialMergedToCertain_array)
    # product after 1 dist merging
    OneDistMergingDf = DataFrameFromOneDistMerging(Certainlist,UncertainReadings,OneCopyDfNotInvolved)
    arrayofInterest = np.concatenate([CertainlistID, UncertainReadingsID, OneCopyReadingsNotInvolved], axis=0)
    # repeat unique and split merging procedure with dist2 and dist3 merging. until line 810.
    TwoDist_pd = CoordsOfTdistInDistMatrix(Mdis,2,arrayofInterest)
    droplist_2distUnique_array,OneDistMergingDf = FindUniqueTwoDistReading(TwoDist_pd,OneDistMergingDf)
    Dist2UncertainReadings,Dist2UncertainReadingsID = DefineUncertainReadings(OneDistMergingDf,5,10)
    Dist2Certainlist,Dist2CertainlistID = DefineCertainlist(OneDistMergingDf,10)
    d2UncertainReadings,d2droplist_PartiallyMergedToUncertain_array = d2Split_Merge_1copy_to_uncertain(Dist2UncertainReadings,Dist2UncertainReadingsID,OneCopyReadingsNotInvolved,Mdis,OneDistMergingDf,arrayofInterest,2)
    d2Certainlist,d2droplist_PartialMergedToCertain_array= d2Split_Merge_1copy_to_certain(Dist2Certainlist,Dist2CertainlistID,OneCopyReadingsNotInvolved,Mdis,OneDistMergingDf,arrayofInterest,2)
    d2OneCopyDfNotInvolved,d2OneCopyReadingsNotInvolved =  OneCopyreadingNotJoinedMerging(OneCopyDfNotInvolved,OneCopyReadingsNotInvolved,d2droplist_PartiallyMergedToUncertain_array, d2droplist_PartialMergedToCertain_array)
    TwoDistMergingDf = DataFrameFromOneDistMerging(Dist2Certainlist,Dist2UncertainReadings,d2OneCopyDfNotInvolved)
    arrayofInterest = np.concatenate([Dist2CertainlistID, Dist2UncertainReadingsID, d2OneCopyReadingsNotInvolved], axis=0)
    ThreeDist_pd = CoordsOfTdistInDistMatrix(Mdis,3,arrayofInterest)
    droplist_3distUnique_array,ThreeDistMergingDf = FindUniqueTwoDistReading(ThreeDist_pd,TwoDistMergingDf)
    Dist3UncertainReadings,Dist3UncertainReadingsID = DefineUncertainReadings(TwoDistMergingDf,5,10)
    Dist3Certainlist,Dist3CertainlistID = DefineCertainlist(TwoDistMergingDf,10)
    droplist_3distUnique_array,ThreeDistMergingDf = FindUniqueTwoDistReading(ThreeDist_pd,TwoDistMergingDf)
    d3UncertainReadings,d3droplist_PartiallyMergedToUncertain_array = d2Split_Merge_1copy_to_uncertain(Dist3UncertainReadings,Dist3UncertainReadingsID,d2OneCopyReadingsNotInvolved,Mdis,TwoDistMergingDf,arrayofInterest,3)
    d3Certainlist,d3droplist_PartialMergedToCertain_array= d2Split_Merge_1copy_to_certain(Dist3Certainlist,Dist3CertainlistID,d2OneCopyReadingsNotInvolved,Mdis,TwoDistMergingDf,arrayofInterest,3)
    d3OneCopyDfNotInvolved,d3OneCopyReadingsNotInvolved =  OneCopyreadingNotJoinedMerging(d2OneCopyDfNotInvolved,d2OneCopyReadingsNotInvolved,d3droplist_PartiallyMergedToUncertain_array, d3droplist_PartialMergedToCertain_array)
    # prepare dataframe of uncertain+certains readings(array2) and ceratain readings(array3)
    seqIDarray2 = np.concatenate([np.array(list(d3Certainlist['SeqID'])),np.array(list(d3UncertainReadings['SeqID']))],axis =0)
    CNumarray2 = np.concatenate([np.array(list(d3Certainlist['C_Num'])),np.array(list(d3UncertainReadings['C_Num']))],axis =0)
    seqIDarray3 = np.concatenate([np.array(list(d3Certainlist['SeqID']))],axis =0)
    CNumarray3 = np.concatenate([np.array(list(d3Certainlist['C_Num']))],axis =0)
    # Functions to get N masked revsed seq with switch ###
    ###array_rev = MaskNwithXDist_pd(OneDist_pd,array_rev,seqIDarray2)
    seq_array = StringMatrixToSeq(array_rev)
    seq_df = pd.DataFrame(seq_array)
    MergedSeqRecords = seq_df.iloc[seqIDarray2]
    MergedSeqRecords_array = np.array(list(MergedSeqRecords[0]))
    MergedSeqRecords2 = seq_df.iloc[seqIDarray3]
    MergedSeqRecords_array2 = np.array(list(MergedSeqRecords2[0]))
    FastaRecords_Final_CandUnC(MergedSeqRecords_array)
    FastaRecords_Final_Certain(MergedSeqRecords_array2)
    ExportMergeFinal(d3Certainlist,ThreeDistMergingDf)


    stop = timeit.default_timer()
    execution_time = stop - start

    print("Program Executed in "+str(execution_time))

   
if __name__ == "__main__":
    main()
