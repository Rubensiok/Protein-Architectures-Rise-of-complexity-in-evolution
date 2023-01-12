#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:55:53 2023

@author: ruben
"""
#%%
import scipy
from tqdm import tqdm
import re
from functools import partial
from nltk import ngrams
from nltk.tokenize import word_tokenize
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer, HashingVectorizer
import plotly.figure_factory as ff
import scipy.cluster.hierarchy as shc
from scipy.cluster import hierarchy
import pandas as pd
from tqdm import tqdm
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, PhyloTree, CircleFace, TextFace,layouts, AttrFace,faces
import pickle
import matplotlib
import random
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import numpy as np

#infile = '/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/todas_arq_comb_pruebaaassssssss.csv'

def get_example_tree(tr):
        # Create a random tree and add to each leaf a random set of motifs
        # from the original set
    t = Tree(tr,format=1)

    for leaf in t.iter_leaves():
        leaf.add_feature('phylum', True)
        # leaf.add_feature('count',True)

    # lis_leaf = []
    for leaf in t.iter_leaves():
        # lis_leaf.append(leaf.name)
        leaf.phylum = dic_phyl[leaf.name]
        # leaf.count = 1
    # print(t.get_ascii(show_internal=True))

    def name_internal_nodes(T):
        for node in T.traverse():
            node.add_feature('count',0)
            if node.is_leaf()==False:
                
                # node.add_feature(count=,True)
                # list names of leaves
                leaf_phyl=[leaf.phylum for leaf in node.iter_leaves()]
                names = list(leaf_phyl)
                names_unique = list(set(leaf_phyl))
                # if all leaves have the same name, give that name to the node
                if (len(names_unique)==1):
                    node.name = names_unique[0]  
                    
                    node.count += len(list(leaf_phyl))
            else:
                node.count = 0
                    # print(node.name,node.count)
    name_internal_nodes(t)

    def collapsed_leaf(node):

        if len(node2labels[node]) == 1:

            return True
        else:
            return False
    node2labels = t.get_cached_content(store_attr=["phylum",'count'])
    nw = t.write(features=['count'],is_leaf_fn=collapsed_leaf)
    t2 = Tree(nw)
    # node.name
    lis_quit = []
    for node in t2.traverse():
        if node.is_leaf()==True:
            if hasattr(node,"count"):
                try:
                    if node.name.replace('"','') in dic_col:
                        node.img_style['fgcolor'] = dic_col[node.name.replace('"','')] 
                        node.img_style['size'] = 5 + (int(node.count)*20/63)
                    elif dic_phyl[node.name.replace('"','')] in dic_col:
                        node.img_style['fgcolor'] = dic_col[dic_phyl[node.name.replace('"','')]]
                        node.img_style['size'] = 5
                        lis_quit.append(node.name)
                except:
                    node.img_style['size'] = 5
            if node.name in grac:
                node.img_style['bgcolor'] = "LightSteelBlue"
            elif node.name in terra:
                node.img_style['bgcolor'] = "DarkSeaGreen"
            else:
                continue
    return t2, lis_quit

def get_example_tree_nocolapse():
        # Create a random tree and add to each leaf a random set of motifs
        # from the original set
    t = Tree(finaltree,format=1)

    for leaf in t.iter_leaves():
        leaf.add_feature('phylum', True)

    # lis_leaf = []
    for leaf in t.iter_leaves():
        # lis_leaf.append(leaf.name)
        leaf.phylum = dic_phyl[leaf.name]
    for leaf in t:
        try:
            leaf.img_style['fgcolor'] = dic_col[leaf.name] 
            leaf.img_style['size'] = 10
        except:
            leaf.img_style['fgcolor'] = dic_col[dic_phyl[leaf.name]]
            leaf.img_style['size'] = 5
            lis_quit.append(leaf.name)
    
    return t


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

#%%
infile = '/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/todas_arq_comb_definitive09112022.csv'

grac = ['Proteobacteria','Spirochaetota','Bacteroidota','Acidobacteriota',
            'Planctomycetota','Verrucomicrobiota','Chlamydiota','Bdellovibrionota','Campylobacterota','Desulfobacterota']
terra = ['Actinobacteriota','Firmicutes','Cyanobacteria','Chloroflexota',
         'Aquificota','Deinococcota','Synergistota','Fusobacteriota','Thermotogota']



org_n_grams = {}
with open (infile,'r') as f:
    # lis_ngrams = []
    for line in tqdm(f):
        # print(line)
        wline = line.rstrip('\n').split(',')
        arq = re.sub('\.[0-9]*','','-'.join(wline[6].split('-')))
        org = wline[2]
        arq_padd = 'NNN-' + arq + '-CCC'
        if org not in org_n_grams:
            org_n_grams[org]=''
            # org_n_grams[org]=[]
        org_n_grams[org] += '-' + arq_padd
        
with open (infile,'r') as f:
    org_king = {}
    # dic_clan = {}
    for line in tqdm(f):
        wline = line.rstrip('\n').split(',')
        # print(wline)
        if wline[2] not in org_king:# or wline[2] not in dic_clan:
            org_king[wline[2]] = []
            # dic_comb[wline[2]] = []
        if wline[4] not in org_king[wline[2]]:
            org_king[wline[2]]=wline[4]

           
with open (infile,'r') as f:
    org_phyl = {}
    # dic_clan = {}
    for line in tqdm(f):
        wline = line.rstrip('\n').split(',')
        if wline[2] not in org_phyl:# or wline[2] not in dic_clan:
            org_phyl[wline[2]] = []
            # dic_comb[wline[2]] = []
        if wline[3] not in org_phyl[wline[2]]:
            org_phyl[wline[2]]=wline[3]

df = pd.DataFrame.from_dict(org_n_grams,orient='index')       
df['phylum'] = ''        
df['kingdom'] = ''
for i in df.index:
    df.kingdom[i] = org_king[i]
    df.phylum[i] = org_phyl[i]

with open (infile,'r') as f:
    lis_unigram = []
    for line in tqdm(f):
        wline = line.rstrip('\n').split(',')
        arq = re.sub('\.[0-9]*','','-'.join(wline[6].split('-')))
        for i in arq.split('-'):
            lis_unigram.append(i.lower())
lis_unigram = list(set(lis_unigram))

df.columns = ['Arch', 'Phylum', 'Kingdom']
# remove missing values
df = df.dropna()

# encode target label
le = LabelEncoder()
# df['Kingdom'] = le.fit_transform(df['Kingdom'])
df['Phylum'] = le.fit_transform(df['Phylum'])

cv = CountVectorizer(analyzer = 'word',ngram_range=(1,1),vocabulary=lis_unigram)
X_FINAL = cv.fit_transform(df.Arch)
matrix = pd.DataFrame(X_FINAL.toarray(),columns=cv.get_feature_names(),index = df.index)

dic_col = {}
for i in matrix.columns:
    dic_col[i] = '-'.join(i.upper().split(' '))
matrix.rename(dic_col,axis=1,inplace=True)

matrix[matrix>1] = 1
matrix.insert(0,'superkingdom','')
matrix.insert(1,'phylum','')


for i in matrix.index:
    matrix.superkingdom[i] = org_king[i]
    matrix.phylum[i] = org_phyl[i]

colors = list(matplotlib.colors.cnames.values())
dic_col = {}
    
for i in org_phyl:
    dic_col[org_phyl[i]]=random.choice(colors)
    
matrix.index = matrix.index.str.replace('(','_').str.replace(')','_').str.replace(':','_').str.replace('.','_').str.strip(' ')
   
dic_phyl = {}
for i in org_phyl:
        dic_phyl[i] = org_phyl[i]

for i in org_phyl:
    dic_phyl[i.replace('(','_').replace(')','_').replace(':','_').replace('.','_').strip(' ')] = org_phyl[i]

matrix1 = matrix.iloc[:,2:]

dend = shc.linkage(matrix1, method='weighted',optimal_ordering = True,metric='jaccard')#,metric='euclidean'
tree = hierarchy.to_tree(dend,False)
finaltree = getNewick(tree, "", tree.dist, matrix1.index)

if __name__ == '__main__':
    t,lis_quit=get_example_tree(finaltree)
    ts = TreeStyle()
    ts.mode = 'c'
    ts.tree_width = 1
    ts.scale =  20
    ts.force_topology = True

    t.unroot()
    t.show(tree_style=ts)
    

lis_quit
matrix2 = matrix1.drop(lis_quit)
dend1 = shc.linkage(matrix2, method='weighted',optimal_ordering = False,metric='jaccard')#,metric='euclidean'
tree1 = hierarchy.to_tree(dend1,False)

finaltree1 = getNewick(tree1, "", tree1.dist, matrix2.index)

if __name__ == '__main__':
    t,lis_quit=get_example_tree(finaltree1)
    ts = TreeStyle()
    ts.mode = 'c'
    ts.tree_width = 1
    ts.scale =  20
    ts.force_topology = True
    t.unroot()
    t.show(tree_style=ts)
    
