import fpGrowth
from math import log2


# get the bacterias that belong to 'extreme' and 'host' habitats from 'bactTaxa_Habitat.txt' file
def create_extreme_and_host_list():
    extreme_list = []
    host_list = []
    with open('bactTaxa_Habitat.txt', newline='') as bact_taxa_file:
        txt_reader = bact_taxa_file.readlines()
        for row in txt_reader:
            lines = row.split(";")
            row_len = len(lines)
            bact_name = lines[1]
            habitat = lines[row_len - 1].rstrip()
            if habitat == 'Extreme Environments':
                extreme_list.append(bact_name)
            if habitat == 'Host':
                host_list.append(bact_name)
    return extreme_list, host_list


#get the COGs each bacteria is consist of out of 'cog_words_bac.txt' file
def get_transactions(extreme_list, host_list):
    cog_dict = {}
    with open('cog_words_bac.txt', newline='') as cog_words_file:
        words_reader = cog_words_file.readlines() 
        transections_dict = {}
        for row in words_reader:
            lines = row.split(" ")
            new_lines = lines[0].split("#")
            bact_name = new_lines[3]
            words = new_lines[4].split('\t')
            if bact_name in extreme_list or bact_name in host_list:
                if bact_name not in transections_dict:
                    transections_dict[bact_name] = [[],0]
                    if bact_name in extreme_list:
                        transections_dict[bact_name][1] = 1
                    if bact_name in host_list:
                        transections_dict[bact_name][1] = 0
                new_cog = []
                for word in words[1:]:
                    word = word.rstrip()
                    if word != 'X' and word != "" and word != " " and word not in transections_dict[bact_name][0]:
                        transections_dict[bact_name][0].append(word)
                        new_cog.append(word)
                cog_tuple =  tuple(new_cog)
                if cog_tuple in cog_dict:
                    cog_dict[cog_tuple] = cog_dict[cog_tuple] + 1
                else:
                    cog_dict[cog_tuple] = 1
    return transections_dict, cog_dict


# entropy calculation
def entropy(class0, class1):
    if class0 <= 0 or class1 <= 0:
        return 0
    return -(class0 * log2(class0) + class1 * log2(class1))


#calculate num of appearances of label 0 in and label 1 in transactions and their entropy
def calculate_all_DB_entropy(transactions):
    label_0_num = 0
    label_1_num = 0
    for t in transactions:
        if transactions[t][1] == 1:
            label_1_num = label_1_num + 1
        if transactions[t][1] == 0:
            label_0_num = label_0_num + 1
    sum_all_label = label_0_num + label_1_num
    class0 = label_0_num / sum_all_label
    class1 = label_1_num / sum_all_label
    s_entropy = entropy(class0, class1)
    return s_entropy, label_0_num, label_1_num


#IG value calculation for itemset
def calculate_IG_value(itemSet, transactions, s_entropy, label_0_num, label_1_num):
    num_of_transactions = label_0_num + label_1_num

    class0_contain_itemSet = 0
    class1_contain_itemSet = 0

    # 1. counting appearances of itemset in class0 and in class1
    for bact in transactions:
       if itemSet.issubset(transactions[bact][0]):
           if transactions[bact][1] == 1:
               class1_contain_itemSet = class1_contain_itemSet + 1
           elif transactions[bact][1] == 0:
               class0_contain_itemSet = class0_contain_itemSet + 1
    # 2. respectively calculate the non appearances of itemset in class 0 and class 1
    class0_NotContain_itemSet = label_0_num - class0_contain_itemSet
    class1_NotContain_itemSet = label_1_num - class1_contain_itemSet

    #3. entropy calculation for number of appearances of itemset in classes
    s1_all = class0_contain_itemSet + class1_contain_itemSet
    #input validation
    if s1_all == 0:
        return 0
    s1_class0 = class0_contain_itemSet / s1_all
    s1_class1 = class1_contain_itemSet / s1_all
    s1_entropy = entropy(s1_class0, s1_class1)

    #4. entropy calculation for number of non appearances of itemset in classes
    s2_all = class0_NotContain_itemSet + class1_NotContain_itemSet
    s2_class0 = class0_NotContain_itemSet / s2_all
    s2_class1 = class1_NotContain_itemSet / s2_all
    s2_entropy = entropy(s2_class0, s2_class1)

    # 5. put it all together for calculate of the information gain
    gain = s_entropy - (s1_all / num_of_transactions * s1_entropy + s2_all / num_of_transactions * s2_entropy)
    return gain


# creating a list containing of all
def creatDataSet(transactions):
    newDataSet = []
    for t in transactions:
        newDataSet.append(transactions[t][0])
    return newDataSet


def find_distinguishing_itemsets(transactions, max_items_list, minSup , cog_dict):
    # 1. arrange the transactions as proper input for createInitSet function 
    dataSet = creatDataSet(transactions)
    # 2. arrange the data as proper input for fptree creation
    initSet = fpGrowth.createInitSet(dataSet)
    # 3. crate fptree
    myFPtree, myHeaderTab = fpGrowth.createTree(initSet, minSup)
    if myHeaderTab is None: # tree is empty
        return max_items_list
    freqItems = []
    # 4. apply fpgrowth algorithm for finding frequent itemsets in data. Eventually [freqitems] contains them
    fpGrowth.mineTree(myFPtree, myHeaderTab, minSup, set([]), freqItems, cog_dict)
    max_IG = 0
    max_item = ""
    # 5. calculate entropy values for whole data
    s_entropy, label_0_num, label_1_num = calculate_all_DB_entropy(transactions)
    # 6. calculate IG score for each itemset while finding the one with the highest IG score value
    for item in freqItems:
        Ig = calculate_IG_value(item, transactions, s_entropy, label_0_num, label_1_num)
        if Ig >= max_IG:
            max_item = item
            max_IG = Ig

    str_max_item = itemToStr(max_item)
    max_items_list[str_max_item] = [max_IG]
    copy_transactions = transactions.copy()

    # 7. find all teansactions in data that contain [max_item] and remove them from data base
    for bact in copy_transactions:
        if max_item.issubset(transactions[bact][0]):
            label = transactions[bact][1]
            max_items_list[str_max_item].append([bact,label])
            del transactions[bact]

    # 8. repeat until no transactions left in data base
    if len(transactions) > 0:
        find_distinguishing_itemsets(transactions, max_items_list, minSup, cog_dict)
    else:
        return max_items_list


def itemToStr(max_item):
    row = ''
    for i in max_item:
        row = row + i + " "
    return row


def print_final(max_items_list):
    count_0 = 0
    count_1 = 0
    for item in max_items_list:
        print( 'Itemset : '+ item + " IG : " + str(max_items_list[item][0]))
        for i in range(1,len(max_items_list[item])):
            bact = max_items_list[item][i][0]
            habitat = max_items_list[item][i][1]
            if habitat == 0 :
                count_0 +=1
            else:
                count_1 +=1
            print(' Bact : ' + bact + ' habitat : ' + str(habitat))
        print(str(count_0+count_1))
        print("num of bact in label 0 :"+ str(count_0))
        print("num of bact in label 1 :"+ str(count_1))


def main():
    # 1. get relevant bacterias lists
    extreme_list, host_list = create_extreme_and_host_list()
    # 2. get a dictionary with bacteria's name as key, containg a list of cogs each bacteria consists of and its class number
    transactions , cog_dict= get_transactions(extreme_list, host_list)
    # 3. define min support value for fpgrowth algorithm to get only itemsets with frequent>=[min_sup]
    min_sup = 200
    max_items_list = {}
    # 4. algorithm for finding distinguishing_itemsets, eventually can be find in [max_items_list]
    find_distinguishing_itemsets(transactions, max_items_list, min_sup , cog_dict)
    print_final(max_items_list)

if __name__ == '__main__':
    main()





