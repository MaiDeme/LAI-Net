import pandas as pd
import os
import subprocess
import numpy as np
import random

#list of functions for extracting training and query sets from reference files

def naming_file(isquery,isreal,nsamplepop,npop,nSNP,data_folder='./data/input/'):
    
    nsample=int(npop*nsamplepop)

    if isreal :
        prefix=data_folder+'real/'+nSNP+'_SNP/'
    else :
        prefix = data_folder+'artificial/'+nSNP+'_SNP/'
    
    if isquery:
        #query file
        name_file= prefix+str(npop)+'pop_'+str(nsample)+'test.vcf.gz'
        name_list_file= prefix+str(npop)+'pop_'+str(nsample)+'test.tsv'
        name_map_file= prefix+str(npop)+'pop_'+str(nsample)+'test.smap'
 
    else :
        #sample file for training
        name_file=prefix+str(npop)+'pop_'+str(nsample)+'train.vcf.gz'
        name_list_file=prefix+str(npop)+'pop_'+str(nsample)+'train.tsv'
        name_map_file = prefix+str(npop)+'pop_'+str(nsample)+'train.smap'   
    
    return name_file,name_list_file,name_map_file,prefix

def all_target_counts_reached(selected_counts, target_count):
    '''
    Function to check if the target counts have been reached for all categories in both groups
    '''
    return all(selected_counts[pop] == target_count for pop in selected_counts)

def extract_from_same_file(reference_file,reference_map_file,target_count_q,target_count_s,sample_map_file,query_map_file,sample_list_file,query_list_file,sample_file,query_file,listpop='all'):
    '''
    Function to select a list for the query and the reference sets from the same file.
    The same sample cannot be in both file and population are in the same amount.
    '''
    reference_map = pd.read_csv(reference_map_file, sep="\t",header=None)

    if (listpop=='all'):
        listpop=reference_map[1].unique()
   
    #Grouping by population 
    grouped = reference_map.groupby(1)
    
    # Initialize two empty DataFrame to store the samples
    selected_samples = pd.DataFrame()
    selected_query=pd.DataFrame()
    
    # Initialize dictionaries to keep track of selected counts for each population category s=sample for training and q=query
    selected_counts_s = {pop: 0 for pop in listpop}
    selected_counts_q = {pop: 0 for pop in listpop}

    # Create a list to keep track of selected samples
    selected_indices = []
    error=False
    print("Randomly selecting samples...\n")
    # Continue selecting samples until the target counts are reached for all categories in both groups
    while not all_target_counts_reached(selected_counts_s, target_count_s) or not all_target_counts_reached(selected_counts_q, target_count_q):

        # Randomly select an index for a sample that hasn't been selected yet
        available_indices = [index for index in reference_map.loc[reference_map.iloc[:, 1].isin(listpop)].index if index not in selected_indices]
        if not available_indices:
            error=True
            print("There is not enough samples for creating the files\n")
            print(selected_counts_q)
            print(selected_counts_s)
            break  # All available samples have been selected

        index = random.choice(available_indices)

        # Select the sample at the chosen index
        row = reference_map.iloc[index].to_frame().T.reset_index(drop=True)
        pop = row.values[0,1]
        if selected_counts_s[pop] < target_count_s:
            selected_samples = pd.concat([selected_samples,row], ignore_index=True)
            selected_counts_s[pop] += 1
        elif selected_counts_q[pop] < target_count_q:
            selected_query = pd.concat([selected_query,row], ignore_index=True)
            selected_counts_q[pop] += 1

        # Update the list of selected sample indices
        selected_indices.append(index)

    if not error:
        #Saving data in separates files
        columns=["Sample","Population code"]
        selected_samples.columns=columns
        selected_query.columns=columns

        np.savetxt(sample_map_file, selected_samples, delimiter="\t", fmt="%s")
        np.savetxt(query_map_file, selected_query, delimiter="\t", fmt="%s")
        print(selected_counts_q)
        print(selected_counts_s)
        print('\n')
        
        np.savetxt(sample_list_file, selected_samples["Sample"], delimiter="\t", fmt="%s")
        np.savetxt(query_list_file, selected_query["Sample"], delimiter="\t", fmt="%s")
        subset_cmd_q = "bcftools view" + " -S " + sample_list_file + " -Ov -o " + sample_file + " " + reference_file
        subset_cmd_s = "bcftools view" + " -S " + query_list_file + " -Ov -o " + query_file + " " + reference_file

        print("Running in command line: \n\t", subset_cmd_q,"\n\t",subset_cmd_s)
        os.system(subset_cmd_q)
        os.system(subset_cmd_s)
        print("\nDone")
    return


def extract_from_file(reference_file,reference_map_file,target_count,map_file,list_file,file,listpop='all'):
    
    #Extracting the samples
    reference_map = pd.read_csv(reference_map_file, sep="\t",header=None)
    
    
    if (listpop=='all'):
        listpop=reference_map[1].unique()
   
    
    #Grouping by population 
    grouped = reference_map.groupby(1)
    
    #listpop=['AFR','EUR','SAS','EAS','AMR']

    # Initialize one DataFrame to store the samples
    selected = pd.DataFrame()

    # Initialize dictionaries to keep track of selected counts for each population category
    selected_counts = {pop: 0 for pop in listpop}

    # Create a list to keep track of selected samples
    selected_indices = []
    error=False
    # Continue selecting samples until the target counts are reached for all populations
    while not all_target_counts_reached(selected_counts, target_count):

        # Randomly select an index for a sample that hasn't been selected yet
        available_indices = [index for index in reference_map.loc[reference_map.iloc[:, 1].isin(listpop)].index if index not in selected_indices]
        if not available_indices:
            error=True
            print("There is not enough samples for creating the files\n")
            print(selected_counts)

            break  # All available samples have been selected

        index = random.choice(available_indices)

        # Select the sample at the chosen index
        row = reference_map.iloc[index].to_frame().T.reset_index(drop=True)
        pop = row.values[0,1]
        if selected_counts[pop] < target_count:
            selected= pd.concat([selected,row], ignore_index=True)
            selected_counts[pop] += 1

        # Update the list of selected sample indices
        selected_indices.append(index)
        
    if not error:
        #Saving data in separates files
        columns=["Sample","Population code"]
        selected.columns=columns
        np.savetxt(map_file, selected, delimiter="\t", fmt="%s")
        
        print(selected_counts)
        print('\n')
        
        np.savetxt(list_file, selected["Sample"], delimiter="\t", fmt="%s")
        
        subset_cmd= "bcftools view" + " -S " + list_file + " -Ov -o " + file + " " + reference_file
        print("Running in command line: \n\t",subset_cmd)
        os.system(subset_cmd)
        print("Done")

def merge_intra_continent(reference_file,reference_map_file,target_count,map_file,list_file,file,special_sample,n_special_sample=None,listpop='all'):
    
    '''
    Given a file, list of population and counts to reach, return a test file.
    If a special sample (Name of a sample from a specific subpop and number) is given it returns a vcf with the reach number 
    '''
   
    #Extracting the samples
    reference_map = pd.read_csv(reference_map_file, sep="\t",header=None)
    
    if (special_sample!= None):
        
        
        special_sample=pd.read_csv(special_sample,sep='\t',header=None)
        special_sample=special_sample.random(n_special_sample)
        #removing the special sample from the the reference_map
        refence_map=reference_map.loc[~reference_map.iloc[1].isin(special_sample.iloc[1].tolist())]

    if (listpop=='all'):
        listpop=reference_map[1].unique()
   
    #Grouping by population 
    grouped = reference_map.groupby(1)
    
    #listpop=['AFR','EUR','SAS','EAS','AMR']

    # Initialize one DataFrame to store the samples
    selected = pd.DataFrame()

    # Initialize dictionaries to keep track of selected counts for each population category
    selected_counts = {pop: 0 for pop in listpop}

    # Create a list to keep track of selected samples
    selected_indices = []
    error=False
    # Continue selecting samples until the target counts are reached for all populations
    while not all_target_counts_reached(selected_counts, target_count):

        # Randomly select an index for a sample that hasn't been selected yet
        available_indices = [index for index in reference_map.loc[reference_map.iloc[:, 1].isin(listpop)].index if index not in selected_indices]
        if not available_indices:
            error=True
            print("There is not enough samples for creating the files\n")
            print(selected_counts)

            break  # All available samples have been selected

        index = random.choice(available_indices)

        # Select the sample at the chosen index
        row = reference_map.iloc[index].to_frame().T.reset_index(drop=True)
        pop = row.values[0,1]
        if selected_counts[pop] < target_count:
            selected= pd.concat([selected,row], ignore_index=True)
            selected_counts[pop] += 1

        # Update the list of selected sample indices
        selected_indices.append(index)
        
    if not error:
        #Saving data in separates files
        columns=["Sample","Population code"]
        selected.columns=columns
        np.savetxt(map_file, selected, delimiter="\t", fmt="%s")
        
        print(selected_counts)
        print('\n')
        
        np.savetxt(list_file, selected["Sample"], delimiter="\t", fmt="%s")
        
        subset_cmd= "bcftools view" + " -S " + list_file + " -Ov -o " + file + " " + reference_file
        print("Running in command line: \n\t",subset_cmd)
        os.system(subset_cmd)
        print("Done")

    
    
        
def add_result(acc,method,nSNP,n,rpop_train,apop_train,rpop_test,apop_test,file="../data/accuracy.tsv"):
    res=pd.read_csv(file,delimiter='\t',header=0)
    if btrain: 
        btrain='r'
    else :
        btrain ='a'
    if btest: 
        btest='r'
    else: 
        btest='a'
    data_acc=[[method,nSNP,btrain,train,btest,test,acc]]
    
    info=pd.DataFrame(data=data_acc,columns=res.columns)
    res=pd.concat([res,info],axis=0)
    res.to_csv(file,sep='\t',index=False)
    print(res)

def output_name(method,nSNP,ntrain,ntest,pop,train,test,prefix='../data/'):
    if (method =='lainet'):
        prefix=prefix+'output_lainet/'
    else :
        prefix=prefix+'output_gnomix/'
        
    prefix=prefix+f'{nSNP}_SNP/'
    
    prefix=prefix+f'train{ntrain*len(pop)}_'
    if train:
        for i in pop:
            prefix=prefix+'r'+i+'-'
    else:
        for j in pop:
            prefix=prefix+'a'+j+'-'
    prefix=prefix+f'test{ntest*len(pop)}_'
    if test:
        for i in pop:
            prefix=prefix+'r'+i+'-'
    else:
        for i in pop:
            prefix=prefix+'a'+i+'-'
            
    return prefix+'/'

def output_name_mixed(type_f,nSNP,train,real_pop,artificial_pop,ntrain,prefix='../data/input/mixed_set/'):
    '''
    if type_f==0 it's a vcf if it's ==1 it's a tsv and if it's ==2 it's a smap
    '''
    prefix=prefix+f'{nSNP}_SNP/'
    if train:
        base=f'train{ntrain*(len(real_pop)+len(artificial_pop))}_'
        for i in real_pop:
            base=base+'r'+i+'-'
        for j in artificial_pop:
            base=base+'a'+j+'-'
    else:
        base=f'test{ntrain*(len(real_pop))}_'
        for i in real_pop:
            base=base+'r'+i+'-'
    prefix=prefix+base
  
    if (type_f==0):
        prefix=prefix+'.vcf.gz'
    elif(type_f==1):
        prefix=prefix+'.tsv'
    elif(type_f==2):
        prefix=prefix+'.smap'
    
    return prefix,base


