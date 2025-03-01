#Weighted maximum parsimony
#For determining which allele is the ancestral in a phylogenetic tree
import pandas as pd
from collections import Counter, defaultdict

#Import nucleotide table
df = pd.read_csv("ASR_input.csv", sep="\t", na_values=['NA', 'NaN', ''],comment="#")

df_new = pd.DataFrame()
df_new = df.iloc[:,2:] #skip position columns
current_position = df["POS"].values.tolist()
df_rows = df_new.values.tolist()


#nucleotide substitution rate is different: ranking CT > AG > AT,CA,GT,CG
nucleotide_pairs_list = ["CT","TC","AG","GA","AT","TA","CA","AC","GT","TG","CG","GC"]
nucleotide_pairs_cost_list = [0.7, 0.7, .8, .8, 1, 1,1.1,1.1,1.2,1.2,1.3,1.3]
row_counter = -1
final_output = []
complete_output = []

#When interested in the cost for the state of the node above humans, this is the node we start with
for sample in df_rows:
    row_counter = row_counter+1
    final_pair = []
    final_pair_allele = []
    favoured_allele = False

    ###PIPELINE START
    #get possible focal node alleles:
    sample_states = pd.Series(sample) #turn to series so i can do dropna
    sample_states= sample_states.dropna()
    sample_states = set(sample_states)
    for state in sample_states:
        output = []
        modified = False
        for sample_SNP in sample:
            if state == sample_SNP:
                output.append(0)
            elif not isinstance(sample_SNP, str):
                output.append(-1) #maybe change so i dont include this in the sum
            else: #now add costs
                replacement_pair = state+sample_SNP
                #get index of string of the nucleotide substitution
                matching_indices_replacement_pair = nucleotide_pairs_list.index(replacement_pair) if replacement_pair in nucleotide_pairs_list else None
                #get cost frm list
                current_cost = nucleotide_pairs_cost_list[matching_indices_replacement_pair]
                output.append(current_cost)
        #dividing the cost lists
        human_list = output[:8]
        monkey_list = output[8:]

        #HUMAN PART
        human_cost_matrix = [4.6912, 4.6912, 4.7241, 4.8801, 4.8802, 5.0102, 5.1812, 5.2475]
        human_cost_matrix_i = [4.6912, 3.5877, 4.3598, 1, 4.4231, 4.6938, 4.3657]
        
        #PART1: Calculate normal cost without considering merging branches
        human_list_raw = [human_list[i] / human_cost_matrix[i] for i in range(len(human_list))]    
        #remove negative values from missing data
        human_list_raw = [0 if x < 0 else x for x in human_list_raw]
        # print(human_list)
        human_cost = sum(human_list_raw)

        def get_adjacent_duplicates_positions(lst):
            duplicate_groups = []  # List to store lists of positions of adjacent duplicates
            temp_group = []  # Temporary list to hold a group of consecutive duplicates

            for i in range(1, len(lst)):
                if lst[i] > 0 and lst[i] == lst[i - 1]:  # Only consider values > 0 for duplicates
                    if not temp_group:  # Start a new group
                        temp_group.append(i - 1)  # First element of the duplicate pair
                    temp_group.append(i)  # Add the second element (or more)
                else:
                    if temp_group:  # If there is a group of duplicates
                        duplicate_groups.append(temp_group)  # Add the current group to the result
                        temp_group = []  # Reset the temporary group

            # Check if the last group of duplicates was not added
            if temp_group:
                duplicate_groups.append(temp_group)

            return duplicate_groups


        duplicate_positions = get_adjacent_duplicates_positions(human_list[:-2]) #excluding HgA0, because it is directly connected to focal node
        human_list_processing = human_list

        
        #remove the replicate entries' position from the human_cost_matrix
        def remove_positions(L2, l1):
            return [L2[i] for i in range(len(L2)) if i not in l1]

        if len(duplicate_positions)>0:
            for replicates in duplicate_positions:
                intermediate_cost_matrix = remove_positions(human_cost_matrix, replicates)
                human_list_processing = remove_positions(human_list, replicates)
                to_add_humanlist = human_list[replicates[0]] #cost of mut
                if replicates[0]>0:
                    to_add_cost1 = human_cost_matrix_i[replicates[0]-1]
                    intermediate_cost_matrix.insert(replicates[0],to_add_cost1)
                    human_list_processing.insert(replicates[0],to_add_humanlist)

                to_add_cost2 = human_cost_matrix_i[replicates[-1]]
                intermediate_cost_matrix.insert(replicates[0]+1,to_add_cost2)
                human_list_processing.insert(replicates[0]+1,to_add_humanlist)

                
            #Calculate cost
            human_list_processing = [human_list_processing[i] / intermediate_cost_matrix[i] for i in range(len(human_list_processing))]    
            
            #remove negative values from missing data
            human_list_processing = [0 if x < 0 else x for x in human_list_processing]
            human_cost_2 = sum(human_list_processing)
            
            #Consider cost from merging branches and original calculation:
            human_cost_list = [human_cost, human_cost_2]
            human_cost = min(human_cost_list)


        #MONKEY PART
        monkey_cost_matrix = [6.7652,6,6,6.9542,7.1139] 
        #make list with focal node and all alleles observed in the monkeys to test different scenarios
        test_alleles = set(monkey_list)

        monkey_samples= sample[8:]
        monkey_alleles= sample[8:]
        if output[8]==output[9]: #bonobo and chimp equal
            print("both chim + bonobo mut the same")
            monkey_cost_matrix = [6.7652,6.699,6.9542,7.1139] #is one longer than monkeylist, because branch to focal node
            monkey_alleles= sample[9:]
        
        #need to add focal node state in front of monkey_alleles
        monkey_alleles.insert(0,state)
        monkey_samples.append(state) #add focal node state among monkey alleles, in case it's not there (but was only found among human sequences) 
        monkey_samples = pd.Series(monkey_samples)
        root_alleles= monkey_samples.dropna()
        root_alleles = set(root_alleles)
        root_output = []
        for root in root_alleles:
            root_cost = []
            for monkey_SNP in monkey_alleles:
                if monkey_SNP == root:
                    root_cost.append(0)
                elif not isinstance(monkey_SNP, str):
                    root_cost.append(-1)
                else:
                    replacement_pair = root+monkey_SNP
                    matching_indices_replacement_pair = nucleotide_pairs_list.index(replacement_pair) if replacement_pair in nucleotide_pairs_list else None
                    current_cost_m = nucleotide_pairs_cost_list[matching_indices_replacement_pair]
                    root_cost.append(current_cost_m)
            root_cost = [root_cost[i] / monkey_cost_matrix[i] for i in range(len(root_cost))]    
            #remove negative values from list
            root_cost = [0 if x < 0 else x for x in root_cost]
            root_cost_sum = sum(root_cost)
            root_output.append(root_cost_sum)
        monkey_cost = min(root_output)
        cost = human_cost + monkey_cost

        #ADD COSTS
        to_add = str(current_position[row_counter])+","+str(state)+"," + str(cost)
        final_pair.append(cost)
        final_pair_allele.append(state)
        complete_output.append(to_add)
    #pick variant with lowest cost
    index_of_smallest_cost = final_pair.index(min(final_pair))
    lowest_cost_allele = final_pair_allele[index_of_smallest_cost]
    final_output.append(lowest_cost_allele)


#Save the smallest cost alleles per position
smallest_cost_output = pd.DataFrame()
smallest_cost_output["pos"]=current_position
smallest_cost_output["focal_node"]=final_output
smallest_cost_output.to_csv('WMP_Focal_Node_lowcost.csv',sep="\t", index=False)

#also save the entire costs for checking 
all_costs_output = pd.DataFrame()
all_costs_output["pos,allele,data"]=complete_output
all_costs_output.to_csv('WMP_Allcosts_Focal_Node.csv',sep="\t", index=False)


