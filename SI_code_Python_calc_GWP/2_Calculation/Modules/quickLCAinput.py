import brightway2 as bw
import os               
import numpy as np       
import pandas as pd     
import time
from brightway2 import *
from bw_recipe_2016 import add_recipe_2016
import bw2analyzer as bwa
import stats_arrays
import matplotlib.pyplot as plt
import sys

eidb = bw.Database('ecoinvent 3.6 cutoff')

def LCI_todo():
    input_prod = str(input('Enter your input product (ref. product):')) 
    return input_prod


def Sel_prod(eidb):
    all_prod_list = []
    while all_prod_list == []: 
        print("no product selected, or product entered cannot be found from ecoinvent, please enter product:")
        sel_refprod = LCI_todo()
        all_prod_list = [act for act in eidb if act['reference product'] == sel_refprod and "market for" not in act['name']]
        if all_prod_list != []:
            break
    return all_prod_list


#define LCAinput class via manual input "all_prod_list"
class LCAinput_manual:
    all_prod_list = []
    def __init__(self):
        print('Define all_prod_list, starting LCI selection')
    def showLCI(self):
        if self.all_prod_list is not None:
            print(*self.all_prod_list, sep = "\n")
        else:
            print("No product selected, please choose LCI")
    def prodFU(self):
        FU = []
        for prod in self.all_prod_list: 
            FU.append({prod.key:1})
        return FU
    def prodname(self):    
        colname = []
        for prod in self.all_prod_list: 
            colname.append(prod.as_dict()["name"] + "; loc:" + prod.as_dict()["location"] + "; ref.prod:" + prod.as_dict()["reference product"])
        return colname


#only activate when LCAinput() being called
class LCAinput():
    all_prod_list = []
    def __init__(self):
        print("define product by calling .getLCI()")
    def getLCI(self):
        self.all_prod_list = Sel_prod(eidb)
    def showLCI(self):
        if self.all_prod_list is not None:
            print(*self.all_prod_list, sep = "\n")
        else:
            print("No product selected, please choose LCI")
    def prodFU(self):
        FU = []
        for prod in self.all_prod_list: 
            FU.append({prod.key:1})
        return FU
    def prodname(self): 
        colname = []
        for prod in self.all_prod_list: 
            colname.append(prod.as_dict()["name"] + "; loc:" + prod.as_dict()["location"] + "; ref.prod:" + prod.as_dict()["reference product"])
        return colname
    

#final LCA calculations     
def LCAcalc(FU, colname, chosen_methods):
    bw.calculation_setups['sameProdCat'] = {'inv':FU, 'ia':chosen_methods}
    #bw.calculation_setups['sameProdCat']
    ProdCat_final = bw.MultiLCA('sameProdCat')
    df = pd.DataFrame(index=chosen_methods, 
             columns=colname,
             data=ProdCat_final.results.T)
    return(df)
