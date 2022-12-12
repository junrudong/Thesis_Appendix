from Bio import SeqIO
import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from sklearn.metrics import pairwise_distances
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""_This software takes aligned ribosomal gene sequences as input, sequencing error are masked by a cutomized
error rate and partially corrected during merging. The pairwise distance matrix is calculated between readings. The indental sequences are merged followed by merging between highly similar sequences. The outputs
including the pairwise distance after 0 distance merging , the merging outcome(ids and copynumbers) and aligned sequences in fasta file._
Args:
     input_file (filepath): The Fasta file
Output_files:
    1. ceratain readings: pandas csv of ceratain readings:with index,seqID and copy number(C_Num)
    2. high-confidence+tentative readings: pandas csv of ceratain+tentative readings, with index,seqID and copy number(C_Num)
    3. final_merged_high_confidence: fasta: sequences of ceratain readings
    4. final_merged_certanandTentative: fasta: sequences of ceratain and tentative readings
    5. Counter_Pd: pandas Counter csv: all readings with SeqID and C_Num after 0 distance merging
    6. fasta: sequences after 0 distance merging
    7. Mdis: pairwise distance matrix

"""
tree_file = '/data/judong/Python_rRNA_alingmentMergePipeline/Results/InputAndTree/HG001_18S.tree'
def branch_length_distribution(tree_file):
    """
    This function creates a tree. 
    Arges: a phylogenetic tree corresponding to the HG rRNA gene fasta alignment
    return: a phylogenetic tree distance and ETE tree.
    """
    distances = []

    tree = ete3.Tree(tree_file, format=1)

    for node in tree.traverse("postorder"):
        if node.is_root() is not True:
            distances.append(node.dist)
    return distances, tree.


def identify_close_seq(tree, distances):
    """
    This function creates a tree. 
    Args: phylogenetic tree distance and ETE tree
    return: find the phylogenetically related nodes
    """
    count = 0
    variant_nodes =  defaultdict(list)
    interval = st.t.interval(
            alpha=0.95,
            df=len(distances)-1,
            loc=np.mean(distances),
            scale=st.sem(distances))

    low_limit = float(interval[0])
    used = set()

    for node in tree.traverse("preorder"):
        nodes_in_distance = []
        count += 1

        if node.name is "":
            node.name = str(count)

        if node.is_leaf() is not true and node is not used:
            
            leaf_children = node.get_leaves()
            #print(leaf_children)
            dist = 0
            for child in leaf_children:
                distance = tree.get_distance(node, child)
                if distance <= low_limit:
                    nodes_in_distance.append(child)
                else:
                    continue

            #print("leaf",leaf_children)
            #print(nodes_in_distanc
            leaf_sorted = set(leaf_children)
            nodes_sorted = set(nodes_in_distance)
            if len(leaf_sorted)*0.9 > len(nodes_sorted):
                used.add(node)
                anc_node = node
                while anc_node.is_root() is False:

                    anc_node = anc_node.up
                    used.add(anc_node)

            if len(leaf_sorted)*0.9 >= len(nodes_sorted):
                #for x in node.iter_descendants():
                 #   used.add(x)

                for child in nodes_in_distance:

                    variant_nodes[node.name].append(child.name)

                    anc_node = child
                    while anc_node != node:
                        anc_node = anc_node.up
                        used.add(anc_node)

                
    final_parents = variant_nodes.copy()
    for variant in variant_nodes:
        #print(variant, variant_nodes[variant])

        try:
            variant_parent = ((tree&variant).up).name
        except:
            continue

        if variant_parent in variant_nodes:
            del final_parents[variant]

    #print(final_parents)
    return final_parents
def TreeIDlist(final_parents):
    Grouplist = list([
            j
        for i,j in final_parents.items()
    ])
    TreeIDlist  = list([
        i
    for ele in Grouplist
        for i in ele
        ])
    return TreeIDlist
def CreateMatrix():
   """This function creates a matrix from a given fasta file. 

    Args:
        input_file (filepath): The Fasta File
    
    return: M(2d array of stringsy): matrix with M_rows of rows and M_cols of columns,reading ids.
    """



        AcChrList = ['chr13','chr14','chr15','chr21','chr22']
        input_file = '/data/judong/Python_rRNA_alingmentMergePipeline/Results/InputAndTree/HG001_Newbash_18S_rd_rename.fas'
        records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
        ids = list(SeqIO.to_dict(SeqIO.parse(input_file, "fasta")).keys())
        ids = list(
                    ids[i]
            for i in range(0,len(ids))
                if ids[i][0:5] in AcChrList
        )
        array =[
            y
            for i in range(0,len(ids))
                if ids[i][0:5] in AcChrList and ids[i] in TreeIDlist
                        for y in records[ids[i]].seq          
        ]
        array2 =[
            y
            for i in range(0,len(ids))
                if ids[i][0:5] in AcChrList and ids[i] not in TreeIDlist
                        for y in records[ids[i]].seq          
        ]
        m = np.array(array)
        M = np.reshape(m, (-1, len(records[ids[1]].seq)))
        M_rows, M_cols = M.shape
        m2 = np.array(array2)
        M2 = np.reshape(m2, (-1, len(records[ids[1]].seq)))
        M2_rows, M2_cols = M2.shape 
        return M,M_rows,M_cols,M2,M2_rows,M2_cols,ids
def CreateReadingNameArray(ids):
       """create an id array. 

    Args:
        input_file (ids): reading ids
    
    return: M(2d array of stringsy): array of ids.
    """
    ReadingNameArray = np.array([ids[i] for i in range(0,len(ids))])
    return ReadingNameArray
def RemoveTerminalGaps(M,M_rows, M_cols):
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


    #the rows with 1st 3 chars are gaps
    MissInfoRowCoordsF=[
        (i)
    for i,j in gapCoords
        if j==0 and M[i,(j+1)]=='-' and M[i,(j+2)]=='-'
]    

    #the rows with last 3 chars are gaps
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
def CreateReadingNameDicts(ReorderRN):
    """
    This function creates a reading name dictionary. 
    Arges(string_array): reading name array
    return: reading name dictionary.
    """
    ReadingDicts = {}
    for i in range(0,len(ids)):
        ReadingDicts[i] = [ReorderRN[i]]
    return ReadingDicts
def CreateReadingNameDicts(ReorderRN):
    ReadingDicts = {}
    for i in range(0,len(ids)):
        ReadingDicts[i] = [ReorderRN[i]]
    
    return ReadingDicts
def ErrorMasking(M,M_rows,M_cols,T):
        """This function masks the error reads by analyzing the majority of readings on the same positions, error rate is  
    one of the input values.

    Args:
        M (2d array of stringsy):Sequence matrix
        M_rows (int): :rows of array of strings
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

        M2[M2 == 'a']=int("1")
        M2[M2 == 'c']=int("2")
        M2[M2 == 't']=int("3")
        M2[M2 == 'g']=int("4")
        M2[M2 == 'n']=int("0")
        M2[M2 == '-']=int("5")#-2
        M2[M2 == '+']=int("6")#-1
        Int = np.vectorize(lambda x: int(x))
        M2 = Int(M1)   
        M2[M2 == 5]=-2
        M2[M2 == 6]=-1
        print("First three position of number Matrix:",M2[0,0],M2[0,1],M2[0,2],"Datatype:",M2.dtype)
        return M2
def CalculatePairWiseDistance_sk(u,v
    """_customizing calculation formula for sklearn pairwise distance matrixfunction_

    Args:
        u (1d array of strings): 1st row of seq reading
        v (1d array of strings): the row of seq reading to pair with 'u' for distance calculation

    Returns:
        a+b(int):the distanct between u and v
    """   
    a = np.bitwise_and(np.bitwise_and(u > 0,v > 0),(u!=v)).sum(axis=0)
    b = np.bitwise_and(u == -2,v>=0).sum(axis=0) + np.bitwise_and(v == -2,u>=0).sum(axis=0)
    return a+b
def CalculatePairWiseDistance(M2,ReadingD):
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
        # distance matrix with M_rows* M_rows shape
        Mdis0 = np.full((M_rows,M_rows),-1)
        repeat = np.full(M_rows,1,dtype=bool)
        counter = np.full(M_rows,1,dtype=int)
        nmcounter = 0
        #1.1 If we want to see How Merge happens unhash the code below
        #MergeDicts_Pd = pd.DataFrame( {'SeqID_i':[],  'SeqID_j':[]})
        for i in range(0,M_rows-1):
            if repeat[i]==False:
                    continue
            else:
                for j in range(i+1,M_rows):
                    if repeat[j]==False:
                         continue
                    else:
                        a = np.bitwise_and(np.bitwise_and(M1[i,]>0,M1[j,]>0),(M1[i,]!=M1[j,])).sum(axis=0)
                        b = np.bitwise_and(M1[i,]== -2,M1[j,]>=0).sum(axis=0) + np.bitwise_and(M1[j,]== -2,M1[i,]>=0).sum(axis=0)
                        Mdis0[i,j]= a+b  
                        if a+b == 0:
                            counter[i]+=1
                            counter[j]-=1
                            repeat[j]= False
                            for ele in ReadingD[j]:
                                if ele not in ReadingD[i]:
                                    ReadingD[i].append(ele)
                            #1.1 If we want to see how merging happens, unhash the code below
                            #df2 = pd.DataFrame([[i,j]], columns=['SeqID_i','SeqID_j'])
                            #MergeDicts_Pd = pd.concat([df2, MergeDicts_Pd])
                      #  else:
        return counter,Mdis0,ReadingD#MergeDicts_Pd
def n_correction(Mdis0,M1):    
    """ This function correct the terminal gaps and the 'n's in by 0 dist consensus.
    arges(2D_int distance matrix, and 2D_in array): distance matrix0, int matrix
    output: int matrix
    
    """
    M_Rows,M_Cols = M1.shape
    def unique_rows(a):
        a = np.ascontiguousarray(a)
        unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
        return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    
    r,c = np.where(Mdis0 == 0)
    XCoords=np.vstack((r,c)).T
    XCoords = unique_rows(np.sort(XCoords))
    for (x,y) in XCoords:
        for i in range(0,M_Cols-1):
            if M1[x,i]==0 or M1[x,i]==-1: 
                M1[x,i] = M1[y,i]
    return M1
def reduction_Matrix(M2,counter):
    """_This function remove the merged readings fromt he 2d array of the readings in int dtype based on the counter.  
    The readings with positive copynumbers were returned as output_

    Args:
        M2 (_2d array of int_): 
        counter (_1d array of int_): output of the function'CalculatePairWiseDistance(M2)'

    Returns:
        array2(2d array of int): M2 after removing the merged reading.
    """    
    M_rows,M_cols = M2.shape
    dict1={}
    for i,j in the enumerate (counter):
        dict1[i]=j
    array2 = np.zeros((len(dict1.keys()),M_cols))
    for i in range(0,len(dict1.keys())):
        if dict1[i]==0:
            continue
        else:
            array2[i]=M2[list(dict1.keys())[i]]
    return array2
def ReducedReadingD(counter,ReadingD):
    """This function createt the reading name dictionary
    Args(1d int array, dictinary): copy number counter,Reading name dictionary
    Output: Reading name dictionary with the removal of merged.
    
    """

    ReadingD2 = {}
    a = np.array([
        x
    for x,y in enumerate(counter)
        if y != 0
       
    ])
    for i in a:
        ReadingD2[i]= ReadingD[i]
    return ReadingD2
def ReindexAfterReduction(ReadingD2):
    Reading_Pd = {}
    for x,y in enumerate(ReadingD2.keys()):
        Reading_Pd[x]=ReadingD2[y]
    return Reading_Pd
def Export(Mdis):
    """function that exporting pairwise distance matrix into csv file, more csvs such as merging path can be exported if necssary.

    Args:
        Mdis (2d array of int), 
        
    """ 
    from numpy import savetxt
    savetxt('Dist_Matrix_v2.0_reduction_1e2.csv',Mdis,fmt='%i',delimiter=',')
    #savetxt('Counter1e2.csv',counter,fmt='%i',delimiter=',')
    #savetxt('Merge_list.csv',Merge_list,fmt='%i',delimiter=',')#   continue
def NumberMatrixToString (reduced_array):
    """_This function transfers 2d array of int back to 2d array of strings_

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
    ### input:string matrix in ['a','c','t','g'], output is ['actg'] array
    rows,cols = array_rev.shape
    separator = ''
    seq_array = np.array([
    separator.join(array_rev[i,])
        for i in range(0,rows)
            ])
    return seq_array
def RemoveGap(array_rev):
    """_This function remove the gap lanes/

    Args:
        rev_array (1d array of strings)
    """    
    a,b = array_rev.shape
    array_gap = [
        i 
        for i in range(0,b-1)
            if sum(array_rev[:,i] == '-')== a
    ]
    array_rev = np.delete(array_rev,array_gap,1)
    return array_rev
def FastaRecords(seq_array):
    """_This function writes 1d array of strings into SeqRecords of SeqIO_, and expor a fasta file

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
    SeqIO.write(records, 'Merged_1e2' ,"fasta")
    return 
def GenerateCounterDataFrame(counter):
    """ This function generates a counter dataframe for each reading and with seqID and copy number(C_Num).

    Args:
        counter (1d array of int): Counter of copyNumbers

    Returns:
        Counter_Pd(pandas): DataFrame with Columns: CopyNumber:C_Num and SeqIDs
    """   
    ### Input: Counter that of copyNumbers
    ### NonZeroArray: Counter without 0 copyNumbers
    ### Counter_Pd: DataFrame with Columns: CopyNumber:C_Num and SeqIDs    
    Counnter_Dict = {}
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

    """Delete the isolated single readings that are parallel with any other readings.
        If a reading's at least 5 step away from any of the other readings, then is will be removed

    Args:
        Mdis (2d array of pairwise_distances):
        Counter_Pd (dataframe of seq readings):
        T (int): Distance threshold that defines distant
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
  
    return Counter_Pd_removeDist,DistList_CNumblowT
def CoordsOfOnedistInDistMatrix(Mdis):
        """This function generate a dataframe of seqIDs which are 1 step away from each other.

    Args:
        Mdis (2d array of pairwise_distances)
    Returns:
        OneDist_pd(dataframe of 2 cols seqIDs): dataframe represents the connections in columns of seqIDi and SeqIDj.
    """ 
    def unique_rows(a):
        """_A function that takes an array tuples as input and removes the duplicated tuples within the tuple array_. 
        Reason of implementation: distance coefficients of the shows up repeatably from the squared pairwise distance matix. In order to avoid 
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
def FindUniqueOneDistReading(OneDist_pd,Counter_Pd,ReadingD2):
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
            if len(OneDist_pd.loc[OneDist_pd["Seq_ID_j"]==j,"Seq_ID_i"]) == 1:
                for ele in ReadingD2[j]:
                    if ele not in ReadingD2[i]:
                        ReadingD2[i].append(ele)
                droplist2[j]= j
                Counter_Pd.loc[Counter_Pd['SeqID'] == i,['C_Num']] += int(Counter_Pd.loc[Counter_Pd['SeqID']==j,'C_Num'])
                droplist2_array = np.array(list(droplist2.keys()))
    try:
        return droplist2_array,Counter_Pd,ReadingD2
    except:
        return [],Counter_Pd,ReadingD2
def Dropthereadings(Counter_Pd,droplistT1_array):
    Counter_Pd = Counter_Pd.drop(droplistT1_array,axis=0)
    for i in droplistT1_array:
        Reading_Pd.pop(i)
    return Counter_Pd,Reading_Pd       
def OneCopyReadingList(Counter_Pd):
    OneCopyDf = Counter_Pd[Counter_Pd["C_Num"]== 1]
    OneCopyDfSeqID = np.array(list(OneCopyDf["SeqID"]))
    return OneCopyDf,OneCopyDfSeqID
def DefineTentative(Counter_Pd,T1,T2):
    """_From Counter_Pd this function isolates the tentative readings of interested:
        define tentative readings of interested: the readings with C_Num between threshold T1 and T2. From the dataframe after Uniquemerging.

    Args:
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist (pandas datafrom with seqIDs and C_Num): datafrom with seqIDs and C_Num, updated C_Num from merging
        T1 (int): _lower threshold for tentative readings of interested_
        T2 (in): _upper threshold for tentative readings of interested_

    Returns:
        TentativeReadings(dataframe of tentative readings of Interest): Counter_Pd removed ditant readigns
        TentativeReadingsID(1d array of int):(array of seqIDs of Tentativereadings
    """ 
    #Input: T1: lower bond of Tentativereadings, T2: upperBond of TentativeReadings
    Tentativelist = (Counter_Pd[Counter_Pd["C_Num"] >T1]).loc[(
        Counter_Pd[Counter_Pd["C_Num"] >T1])['C_Num']<T2]
    TentativeID = np.array(list(Tentativelist['SeqID']))
    return Tentativelist,TentativeID

def DefineConfidence(Counter_Pd,T1):
    """_Define Confidence readings list of interest. 
    Confidence Reading: readings with C_Num that above T1_

    Args:
        Counter_Pd_Unique_RemoveUniqueMerge_removeDist (pandas datafrom with seqIDs and C_Num): datafrom with seqIDs and C_Num, updated C_Num from merging
        T1 (int): _upper threshold for tentative readings of interested_

    Returns:
       Confidencelist(dataframe)
       ConfidencelistID(1d array of int)
    """   
    #Input: T1 defines the lower limit for Confidence C_Num
    Confidencelist = Counter_Pd[Counter_Pd["C_Num"] >=T1] 
    ConfidenceID = np.array(list(Confidencelist['SeqID']))
    return Confidencelist,ConfidenceID
def Split_Merge_1copy_to_TentativeTentative(Tentativelist,TentativeID,OneCopyDfSeqID,Mdis,Counter_Pd,Reading_Pd,T):
    """_Split merging readings with 1 copies to tentative readings_

    Args:
        TentativeReadings (_pdDataframe_): 
        TentativeReadingsID (_np1darray_): selected id array
        OneCopyDfSeqID (_np1darray_): selected id array
        Mdis (_2d array of int_): distance matrix
        Counter_Pd (_pdDataframe_): dataframe with SeqID and C_Num
        T (int): distance of of selection

    Returns:
       TentativeReadings(dataframe)
       droplist_PartiallyMergedToTentative_array(1d array of int)
    """    
    Mdis_Pd = pd.DataFrame(Mdis)
    droplist_PartialTentative = {}
    
    for i in TentativeID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                droplist_PartialTentative[j] = j
                droplist_PartialTentative_array = np.array(list(droplist_PartialTentative.keys()))
                Tentativelist.loc[Tentativelist['SeqID'] == i,'C_Num'] +=  float(Tentativelist.loc[Tentativelist['SeqID'] == i,'C_Num']) /sum([
                float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                ])
                for ele in Reading_Pd[j]:
                    if ele not in Reading_Pd[i]:
                           Reading_Pd[i].append(ele)   
    try:    
        return Tentativelist,droplist_PartialTentative_array,Reading_Pd
    except:
        return Tentativelist,[],Reading_Pd
def Split_Merge_1copy_to_confidence(Confidencelist,ConfidenceID,OneCopyDfSeqID,Mdis,Counter_Pd,Reading_Pd,T):
    droplist_PartialConfidence = {}
    Mdis_Pd = pd.DataFrame(Mdis) 
    for i in ConfidenceID:
        for j in list((Mdis_Pd.loc[Mdis_Pd[i]==T]).index):
            if j in OneCopyDfSeqID:
                droplist_PartialConfidence[j] = j
                droplist_PartialConfidence_array = np.array(list(droplist_PartialConfidence.keys()))
                Confidencelist.loc[Confidencelist['SeqID'] == i,'C_Num'] += float(Confidencelist.loc[Confidencelist['SeqID'] == i,'C_Num']) /sum([
                float(Counter_Pd.loc[Counter_Pd['SeqID']==k,'C_Num'])
                for k in list(Mdis_Pd.loc[Mdis_Pd[j]==T].index)
                ])
                for ele in Reading_Pd[j]:
                    if ele not in Reading_Pd[i]:
                        Reading_Pd[i].append(ele)                
    try: 
        return Confidencelist,droplist_PartialConfidence_array,Reading_Pd
    except:
        return Confidencelist,[],Reading_Pd


def OneCopyreadingNotJoinedMerging(OneCopyDf,OneCopyDfSeqID,droplist_PartialTentative_array, droplist_PartialConfidence_array):
"""_Output a dataframe and an array list of 1 copy readings that have not joined any 1 dist merging_

    Args:
        OneCopyDf (dataframe): 
        OneCopyDfSeqID (1d array of seqIDs): 
        droplist_PartiallyMergedToTentative_array (1d array of seqIDs): 
        droplist_PartialMergedToConfidence_array (1d array of seqIDs): 

    Returns:
       OneCopyDfNotInvolved(dataframe)
       OneCopyReadingsNotInvolved(1d array of seqIDs)
    """    
    a = np.concatenate((droplist_PartialTentative_array, droplist_PartialConfidence_array), axis=None)
    ReadingsInvolved = np.unique(a)
    OneCopyNotInvolved_array = OneCopyDfSeqID[~np.isin(OneCopyDfSeqID, ReadingsInvolved)]
    OneCopyNotInvolved_df = OneCopyDf.drop(ReadingsInvolved)
    return OneCopyNotInvolved_df,OneCopyNotInvolved_array
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
def autodistmerge(Mdis,arrayofInterest,OneCopyNotInvolved_df,OneCopyNotInvolved_array,Reading_Pd,ids,Counter
    
    T = 1
    CheckPoint = 1
    while CheckPoint > 0.1:
        T +=1
        TDist_pd = CoordsOfTdistInDistMatrix(Mdis,T,arrayofInterest)
        droplist_Unique_array,Counter_Pd,Reading_Pd = FindUniqueTDistReading(TDist_pd,Counter_Pd,OneCopyNotInvolved_array,Reading_Pd)
        Tentativelist,TentativeID = DefineTentative(Counter_Pd,5,10)
        Confidencelist,ConfidenceID = DefineConfidence(Counter_Pd,10)
        Tentativelist,droplist_PartialTentative_array,Reading_Pd = Split_Merge_1copy_to_Tentative(Tentativelist,TentativeID,OneCopyNotInvolved_array,Mdis,Counter_Pd,Reading_Pd,T)
        Confidencelist,droplist_PartialConfidence_array,Reading_Pd= Split_Merge_1copy_to_Confidence(Confidencelist,ConfidenceID,OneCopyNotInvolved_array,Mdis,Counter_Pd,Reading_Pd,T)
        OneCopyNotInvolved_df,OneCopyNotInvolved_array =  OneCopyreadingNotJoinedMerging(OneCopyNotInvolved_df,OneCopyNotInvolved_array,droplist_PartialTentative_array, droplist_PartialConfidence_array)
        arrayofInterest = np.concatenate([ConfidenceID, TentativeID, OneCopyNotInvolved_array])
        CheckPoint = len(OneCopyNotInvolved_array)/len(ids) 
    else:
        print(len(OneCopyNotInvolved_array),'of One Copy Readinds has not involved in merging,','Merging happened up to',T,'dist')
def FindUniqueTDistReading(T_pd,Counter_Pd,OneCopyNotInvolved_array,Reading_Pd):
    droplist_TUnique = {}
    for i in range(0,len(Counter_Pd)):
        for j in T_pd.loc[T_pd["Seq_ID_i"]==i,"Seq_ID_j"]:
            if len(T_pd.loc[T_pd["Seq_ID_i"]==j,"Seq_ID_j"]) == 1:
                if j in OneCopyNotInvolved_array:
                    droplist_TUnique[j]= j
                    Counter_Pd.loc[Counter_Pd['SeqID'] == i,['C_Num']] += float(Counter_Pd.loc[Counter_Pd['SeqID']==j,'C_Num'])
                    droplist_TUnique_array = np.array(list(droplist_TUnique.keys()))
                    for ele in Reading_Pd[j]:
                        if ele not in Reading_Pd[i]:
                            Reading_Pd[i].append(ele)                
    try:
        return droplist_TUnique_array,Counter_Pd,Reading_Pd
    except:
        return [],Counter_Pd,Reading_Pd
def CalChrCount(ids,Df):    

    UniChrPos = np.array( ['chr13','chr14','chr15','chr21','chr22'])
    IDarray = np.array(Df['SeqID'])
    data = np.full((len(IDarray),len(UniChrPos)),0)
    UniChrPos_df = pd.DataFrame(data,columns=UniChrPos)
    UniChrPos_df['SeqID']=IDarray
    for key in Reading_Pd.keys():
        for x,y in enumerate(Reading_Pd[key]):
            try:
                UniChrPos_df.loc[UniChrPos_df['SeqID']==key,[Reading_Pd[key][x][0:5]]] +=1
            except:
                continue
    IDlist = list(UniChrPos_df['SeqID'])
    UniChrPos_df = UniChrPos_df.set_index('SeqID',drop=True)   
    return IDlist,UniChrPos_df
def ExpChrCount(CU_count,C_count):
    columns_CU = list(C_count.columns)
    columns_C =  list(C_count.columns)
    data_CU = np.full((1,len(columns_CU)),0)
    data_C = np.full((1,len(columns_C)),0)
    EmptyRef_CU = pd.DataFrame(data_CU,columns=columns_CU,index=[111111])
    EmptyRef_C = pd.DataFrame(data_C,columns=columns_C,index=[111111])
    ConfidenceAndUnEx = CU_count.append(EmptyRef_CU)
    ConfidenceEx = C_count.append(EmptyRef_C)
    ConfidenceAndUnEx.to_csv(F"LociPosof{len(ConfidenceAndUnEx)}Readings.csv",sep = ',')
    ConfidenceEx.to_csv(F"LociPosof{len(ConfidenceEx)}Readings.csv",sep = ',')
def FastaRecords_Final_CandUnC(MergedSeqRecords_array,IDlist,S):
    rows, = MergedSeqRecords_array.shape
    input_file = '/data/judong/Bash_RibosomalRnaGeneAlignment/AlignmentResults/Ref/ref.fas'
    records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    ids = list(SeqIO.to_dict(SeqIO.parse(input_file, "fasta")).keys())
    array =[
            y
            for i in range(0,len(ids))
                for y in records[ids[i]].seq
        ]
    rows, = MergedSeqRecords_array.shape
    separator = ''
    a = np.array([MergedSeqRecords_array[i].item() for i in range(0,rows)])
    b = [Seq(i) for i in a]
    c = np.array((IDlist,b),dtype=object).T
    d = np.array([['ref',(records[ids[0]].seq)]],dtype=object)
    e = np.array([['ref',(records[ids[1]].seq)]],dtype=object)
    if S == 18:
        f = np.append(c,d,axis = 0)
    if S == 28:
        f = np.append(c,e,axis = 0)

    
    records2 = [
           SeqRecord(Seq(seq), id = str(index), description = "") 
           for index,seq in f 
    ]
    
    SeqIO.write(records2, 'Final_Merged_ConfidenceandTentative' ,"fasta")
    return 
def ConfidenceC_Numlist(Confidence_Df,ConfidenceIDlist):
    ConfidenceC_Numlist = np.array([
        Confidence_Df['C_Num'][i]
             for i in ConfidenceIDlist   
    ])
    ConfidenceID_CNum = np.array([
    f"ID{ConfidenceIDlist[i]}_C{int(ConfidenceC_Numlist[i])}"
        for i in range(0,len(ConfidenceIDlist))
    ])
    return ConfidenceID_CNum
def FastaRecords_Final_Confidence(Confidence,ConfidenceID_CNum,S):
    """_This function is able to out a Fasta file from cerain list of data of interested _

    Args:
        MergedSeqRecords_array (arrayofMegedseq): 
    """    
    input_file = '/data/judong/Bash_RibosomalRnaGeneAlignment/AlignmentResults/Ref/ref.fas'
    records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    ids = list(SeqIO.to_dict(SeqIO.parse(input_file, "fasta")).keys())
    array =[
            y
            for i in range(0,len(ids))
                for y in records[ids[i]].seq
        ]
    rows, = Confidence.shape
    separator = ''
    a = np.array([Confidence[i].item() for i in range(0,rows)])
    b = [Seq(i) for i in a]
    c = np.array((ConfidenceID_CNum,b),dtype=object).T
    d = np.array([['ref',(records[ids[0]].seq)]],dtype=object)
    e = np.array([['ref',(records[ids[1]].seq)]],dtype=object)
    if S == 18:
        d = np.append(c,d,axis = 0)
    if S == 28:
        d = np.append(c,e,axis = 0)

    
    records2 = [
           SeqRecord(Seq(seq), id = str(index), description = "") 
           for index,seq in d 
    ]
    
    SeqIO.write(records2, 'Final_Merged_Confidence' ,"fasta")
    return 

def ExportMergeFinal(Confidencedf,ThreeDistMergingDf):
    """_Export final dataframes_

    Args:
        Confidencedf (dataframe): 
        ThreeDistMergingDf (dataframe): 
    """   
    Confidencedf.to_csv('ConfidenceReadings.csv',sep=',')
    a = ThreeDistMergingDf.loc[ThreeDistMergingDf['C_Num']>=5,]
    a.to_csv('ConfidenceandTentativeReadings.csv',sep=',')
     

