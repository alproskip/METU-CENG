import sys
import pickle
from graphviz import Digraph
from math import *
import numpy as np

def ent(x, y):
    if x == 0 or y == 0: 
        return 0
    return -1 * (x/(x+y) * log(x/(x+y), 2) + y/(x+y) * log(y/(x+y), 2))

def divide(data, attr_index, attr_vals_list):
    """Divides the data into buckets according to the selected attr_index.
    :param data: Current data in the node
    :param attr_index: Selected attribute index to partition the data
    :param attr_vals_list: List of values that attributes may take
    :return: A list that includes K data lists in it where K is the number
     of values that the attribute with attr_index can take
    """
    values = attr_vals_list[attr_index]
    v_num = len(values)
    buckets = []
    for i in range(v_num):
        buckets.append([x for x in data if x[attr_index]==values[i]])
    return buckets

def entropy(data, attr_vals_list):
    """
    Calculates the entropy in the current data.
    :param data: Current data in the node
    :param attr_vals_list: List of values that attributes may take
    (Last attribute is for the labels)
    :return: Calculated entropy (float)
    """
    vals = attr_vals_list[-1]
    c1 = vals[0]
    c2 = vals[1]
    data_count = len(data)
    acc_count = len([x for x in data if x[-1]==c1])
    unacc_count = len([x for x in data if x[-1]==c2])
    return ent(acc_count, unacc_count)

def info_gain(data, attr_index, attr_vals_list):
    """
    Calculates the information gain on the current data when the attribute with attr_index is selected.
    :param data: Current data in the node
    :param attr_index: Selected attribute index to partition the data
    :param attr_vals_list: List of values that attributes may take
    :return: information gain (float), buckets (the list returned from divide)
    """
    ent_data = entropy(data, attr_vals_list)
    childs = divide(data, attr_index, attr_vals_list)
    data_size = len(data)
    child_sizes = [len(x) for x in childs]

    choices = attr_vals_list[attr_index]
    child_ents = [entropy(c, attr_vals_list) for c in childs]

    ent_attr = 0    
    for i in range(len(childs)):
        ent_attr += (child_sizes[i]/data_size) * child_ents[i]
    gain = ent_data - ent_attr
    return gain, childs

def gain_ratio(data, attr_index, attr_vals_list):
    """
    Calculates the gain ratio on the current data when the attribute with attr_index is selected.
    :param data: Current data in the node
    :param attr_index: Selected attribute index to partition the data
    :param attr_vals_list: List of values that attributes may take
    :return: gain_ratio (float), buckets (the list returned from divide)
    """
    data_size = len(data)
    gain, childs = info_gain(data, attr_index, attr_vals_list)
    child_sizes = [len(x) for x in childs]
    multipliers = [(-x/data_size)*log((x/data_size),2) for x in child_sizes if x > 0]
    intrinsic = sum(multipliers)
    return gain/intrinsic, childs

def gini(data, attr_vals_list):
    """
    Calculates the gini index in the current data.
    :param data: Current data in the node
    :param attr_vals_list: List of values that attributes may take
    (Last attribute is for the labels)
    :return: Calculated gini index (float)
    """
    vals = attr_vals_list[-1]
    c1 = vals[0]
    c2 = vals[1]
    c_size = len(data)
    acc_count = len([x for x in data if x[-1]==c1])
    unacc_count = len([x for x in data if x[-1]==c2])
    return (1 - (acc_count/c_size)**2 - (unacc_count/c_size)**2)

def avg_gini_index(data, attr_index, attr_vals_list):
    """
    Calculates the average gini index on the current data when the attribute with attr_index is selected.
    :param data: Current data in the node
    :param attr_index: Selected attribute index to partition the data
    :param attr_vals_list: List of values that attributes may take
    :return: average gini index (float), buckets (the list returned from divide)
    """
    result = 0
    data_size = len(data)
    childs = divide(data, attr_index, attr_vals_list)
    for child in childs:
        c_size = len(child)
        if c_size == 0: continue
        acc_count = len([x for x in child if x[-1]=='acc'])
        unacc_count = len([x for x in child if x[-1]=='unacc'])
        result += (len(child)/data_size) * gini(child, attr_vals_list)
    return result, childs

def chi_squared_test(data, attr_index, attr_vals_list):
    """
    Calculated chi squared and degree of freedom between the selected attribute and the class attribute
    :param data: Current data in the node
    :param attr_index: Selected attribute index to partition the data
    :param attr_vals_list: List of values that attributes may take
    :return: chi squared value (float), degree of freedom (int)
    """
    classes = attr_vals_list[-1]
    vals = attr_vals_list[attr_index]
    c1, c2 = classes[0], classes[1]
    size = len(data)
    acc_count, unacc_count = len([x for x in data if x[-1]==c1]), len([x for x in data if x[-1]==c2])

    acc_ratio, unacc_ratio = acc_count/size, unacc_count/size
    obs, exp = np.zeros((len(vals),2)), np.zeros((len(vals),2))
    splits = divide(data, attr_index, attr_vals_list)

    for i in range(len(splits)):
        obs[i][0] = len([x for x in splits[i] if x[-1]==c1])
        obs[i][1] = len([x for x in splits[i] if x[-1]==c2])

    for i in range(len(splits)):
        ssize = sum(obs[i])
        exp[i][0] = acc_ratio*ssize
        exp[i][1] = unacc_ratio*ssize

    for i in range(2):
        if obs[:,i].any() == False:
            obs = np.delete(obs, i, 1)
            exp = np.delete(obs, i, 1)
            break
    del_idx = []
    for i in range(len(vals)):
        if obs[i].any() == False:
            del_idx.append(i)
    if len(del_idx)!=0:
        obs = np.delete(obs, del_idx, 0)
        exp = np.delete(exp, del_idx, 0)

    m = (obs-exp)**2/exp
    c = sum(sum(m))
    return c, (obs.shape[1]-1)*(obs.shape[0]-1)


def ID3(data, target_Attribute, attr_names, attr_vals_list, heuristic):
    if len([x for x in data if x[-1]=='unacc']) == 0:
        return ['acc','leaf',(0,len(data))]
    if len([x for x in data if x[-1]=='acc']) == 0: 
        return ['unacc', 'leaf', (len(data),0)]
    if len(attr_names)==0:
        acc_count = len([x for x in data if x[-1]=='acc'])
        if acc_count > len(data)/2:
            return ['acc','leaf', (len(data)-acc_count, acc_count)]
        else:
            return ['unacc', 'leaf', (len(data)-acc_count, acc_count)]
    #BEGIN
    if heuristic != 'gini':
        max_gain = -1
        for a in attr_names:
            if(heuristic=='info'):
                gain, div = info_gain(data, a[0], attr_vals_list)
            elif(heuristic=='gain' or heuristic=='gain-c'):
                gain, div = gain_ratio(data, a[0], attr_vals_list)
            if gain > max_gain:
                max_gain=gain
                divisions = div
                best_attr = a[0]
    else:
        min_gini_idx = 99
        for a in attr_names:
            gini, div = avg_gini_index(data, a[0], attr_vals_list)
            if gini < min_gini_idx:
                min_gini_idx = gini
                divisions = div
                best_attr = a[0]
    if heuristic=='gain-c':
        alphas = [0, 2.71, 4.61, 6.25, 7.78, 9.24, 10.5, 12, 13.4, 14.7, 16]
        chisq, df = chi_squared_test(data, best_attr, attr_vals_list)
        if chisq < alphas[df]:
            acc_count = len([x for x in data if x[-1]=='acc'])
            if acc_count > len(data)/2:
                return ['acc','leaf', (len(data)-acc_count, acc_count)]
            else:
                return ['unacc', 'leaf', (len(data)-acc_count, acc_count)]

    attr = [x for x in attr_names if x[0]==best_attr][0]
    node = [attr]
    for i in attr_names:
        if i[0]==best_attr:
            attr_names.remove(i)

    if attr[1] == 'buying': #['vhigh', 'high', 'med', 'low']
        node = [attr, ['VHIGH', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['HIGH', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['MED', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['LOW', ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node += [[ ['VHIGH'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['HIGH']+[ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['MED']+[ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['LOW']+[ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]]
    elif attr[1] == 'maint': #['vhigh', 'high', 'med', 'low']
        node = [attr, ['VHIGH', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['HIGH', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['MED', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['LOW', ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node += [[ ['VHIGH'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)],['HIGH'] + [ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['MED'] + [ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['LOW'] + [ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]]
    elif attr[1] == 'doors': #['2', '3', '4', '5more']
        node = [attr, ['2', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['3', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['4', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['5MORE', ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node  + = [[ ['2'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)],['3'] + [ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['4'] + [ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)], ['5MORE'] + [ID3(divisions[3], 0, attr_names[:], attr_vals_list, heuristic)]]
    elif attr[1] == 'persons': #['2', '4', 'more']
        node = [attr, ['2', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['4', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['MORE', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node  + = [[ ['2'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)],['4'] + [ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['MORE'] + [ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]]
    elif attr[1] == 'lug_boot': #['small', 'med', 'big']
        node = [attr, ['SMALL', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['MED', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['BIG', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node  + = [[ ['SMALL'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)],['MED'] + [ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['BIG'] + [ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]]
    elif attr[1] == 'safety': #['low', 'med', 'high']
        node = [attr, ['LOW', ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)], ['MED', ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)], ['HIGH', ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]] 
        #node  + = [[ ['LOW'] + ID3(divisions[0], 0, attr_names[:], attr_vals_list, heuristic)],['MED'] + [ID3(divisions[1], 0, attr_names[:], attr_vals_list, heuristic)],['HIGH'] + [ID3(divisions[2], 0, attr_names[:], attr_vals_list, heuristic)]]

    return node

def test(tree, test_data, attr_vals_list):
    acc = 0
    for example in test_data:
        temp_tree = tree[:]
        while True:
            if temp_tree[1]=='leaf':
                result = temp_tree[0]
                break
            attr_index = temp_tree[0][0]
            values = attr_vals_list[attr_index]
            new_idx = values.index(example[attr_index]) + 1
            if len(temp_tree[new_idx]) == 1:
                temp_tree = temp_tree[new_idx][0][1]
            else:
                temp_tree = temp_tree[new_idx][1]
        if result == example[-1]:
            acc += 1
    return acc/len(test_data)*100

def print_tree(node, dot, indices, current_idx):
    if node[1]=='leaf':
        return

    for child in node[1:]:
        indices.append(indices[-1]+1)
        child_idx = indices[-1]
        if child[1][1] == 'leaf':
            child_attr_name = child[1][0]+str(child[1][2])
        else:
            child_attr_name = child[1][0][1]
        dot.node(str(child_idx), child_attr_name)
        dot.edge(str(current_idx), str(child_idx), str(child[0]))
        print_tree(child[1], dot, indices, child_idx)


if __name__ == '__main__':
    
    mode = sys.argv[1]
    with open('hw3_data/dt/data.pkl', 'rb') as f:
        train_data, test_data, attr_vals_list, attr_names = pickle.load(f)

    max_gain = -1
    for i in range(5):
        gain, divisions = info_gain(train_data, i, attr_vals_list)
        if gain > max_gain:
            best_attr = i

    tree = ID3(train_data, best_attr, list(enumerate(attr_names)), attr_vals_list, mode)
    test_result = test(tree, test_data, attr_vals_list)
    print("Accuracy is : %"+"%.1f" % test_result)
    
    try:
        plot = sys.argv[2]
        if mode == 'info': tree_name = "Info_Gain"
        elif mode == 'gini':tree_name = "Gini_Index"
        elif mode == 'gain':tree_name = "Gain_Ratio"
        elif mode == 'gain-c': tree_name = "Gain_Chi"
        dot = Digraph()
        dot.node('1',tree[0][1])
        indices = [1]
        print_tree(tree, dot, [1], 1)
        dot.render(tree_name, view=True)
    except Exception as e:
        pass