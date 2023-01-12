#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 14:15:39 2023

@author: ruben
"""


#%%
"""
###############################################################################
##############                  FUNCIONES TFM                    ##############
###############################################################################
"""
from time import time
import os
from pathlib import Path
import pandas as pd
import re
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import csv
import requests
from multiprocessing.pool import ThreadPool
import plotly
import plotly.express as px

#%%
  
def fetch_url(entry):
    path, uri = entry
    if not os.path.exists(path):
        r = requests.get(uri, stream=True)
        if r.status_code == 200:
            with open(path, 'wb') as f:
                for chunk in r:
                    f.write(chunk)
    return path

def protein_comb_new_NEWDATA_repeats(path,obj1,obj2,database,version):
    """
    
    Parameters
    ----------
    path : string
        Folder to create the file with every architectures
    obj1 : string
        First value of pfam_scan output to take into account (hmm_acc).
    obj2 : string
        First value of pfam_scan output to take into account (clan).
    database : string
        folder with every pfam_scan.dom file already separated in different
        folders by taxonomy.
    version : string
        version of the analysis.

    Returns
    -------
    None.

    """
    inicio = time()
    lista_superr = os.listdir(database)
    carpeta = Path(database)
    directorio = str(carpeta).replace('\\','/')
    nombres_proteomas= []
    tamanio_proteomas = []
    superreino = []
    rep_or_n = {}
                  
    f = open (path + '/todas_arq_comb_' + version + '.csv','w')
    df = pd.read_csv('D:/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/uniprot_dataset_junio.csv') # dataframe with every organisms selected from UNIPROT Reference proteome
    df = df[df.SUPERREGNUM!='viruses']
    for superr in lista_superr:
        lista_phylum = os.listdir(directorio + "/" + superr)   
        for phyl in lista_phylum:
            lista_org = os.listdir(directorio + "/" + superr + '/' + phyl) 
            for org in lista_org:
                arch_dom={}
                arch_comb={}
                tamanio = []
                data = open(directorio + "/" + superr + "/" + phyl + "/" + org, "r")
                nombres_proteomas = df['Species Name\n'][df.Proteome_ID==org.replace('.dom','')].iloc[0].strip()
                superreino.append(superr)
                data2 = []
                for linea in data:
                    if linea.startswith("sp") == True or linea.startswith("tr")==True:
                        filaDelCSV = re.split(" +",linea.strip())
                        filaDelData = {
                            "seqId": filaDelCSV[0],
                            "inicio": filaDelCSV[1],
                            "hmm_acc": re.sub('\.[0-9]*','',filaDelCSV[5]),
                            "domfam": filaDelCSV[6],
                            "tipo": filaDelCSV[7],
                            "clan": filaDelCSV[14],
                            "longitud": int(filaDelCSV[10]),
                            "e_val": float(filaDelCSV[12])
                        }
                        data2.append(filaDelData)
                        tamanio.append(filaDelCSV[0])
                tamanio = len(set(list(tamanio)))
                tamanio_proteomas.append(tamanio)
                data.close()
                
                data2 = sorted(data2, key=lambda k: (k['seqId'].lower(), int(k['inicio'])))
                for linea in tqdm(data2):
                    if linea['e_val'] < 0.01:
                        if linea['seqId'] not in arch_comb:
                            arch_comb[linea['seqId']] = []
                            arch_dom[linea['seqId']] = []
                            rep_or_n[linea['seqId']] = []
                        arch_dom[linea['seqId']].append(linea[obj1])
                        rep_or_n[linea['seqId']].append(linea['tipo'])
                        if linea[obj2]!='No_clan':
                            arch_comb[linea['seqId']].append(linea[obj2])
                        if linea[obj2]=='No_clan':
                            arch_comb[linea['seqId']].append(linea[obj1])
                    else:
                        continue

                for key in arch_comb.keys():
                    f.write(key.split('|')[1] + ',' + org.replace('.dom','') + ',' + nombres_proteomas + ',' + phyl + "," + superr + ',' + str(tamanio) + ',') ##### Había una coma de más....
                    f.write('-'.join(arch_comb[key]) + ',' + '-'.join(arch_dom[key]) + ',' + '-'.join(rep_or_n[key]))
                    f.write('\n')                     
    f.close()  
    fin = time()
    print((fin-inicio)/60)


def get_phyla(dic,lis_order):
    dic_phyl = {}
    for i in dic:
        if org_king[i]=='Bacteria':
            phyl = org_phyl[i]
        else:
            phyl = org_king[i]
        if phyl not in dic_phyl:
            dic_phyl[phyl] = []
        for x in dic[i]:
            dic_phyl[phyl].append(x)
    for i in dic_phyl:   
        dic_phyl[i] = list(set(dic_phyl[i]))
    dic_phyl = {k:dic_phyl[k] for k in lis_order}
    return dic_phyl


def get_phyla_TOT(dic):
    dic_phyl = {}
    for i in dic:
        phyl = org_phyl[i]
        if phyl not in dic_phyl:
            dic_phyl[phyl] = []
        for x in dic[i]:
            dic_phyl[phyl].append(x)
    for i in dic_phyl:   
        dic_phyl[i] = list(set(dic_phyl[i]))
    return dic_phyl


def get_phyla_prot(dic):
    dic_phyl = {}
    for i in dic:
        phyl = org_phyl[i]
        if phyl not in dic_phyl:
            dic_phyl[phyl] = []
        for x in dic[i]:
            dic_phyl[phyl].append(x)
    return dic_phyl
    

def matriz_inter(referencia, dic):
    lista_len = []
    organismos = list(dic.keys())
    for org in organismos:
        prueba = len(set(dic[org]) & set(dic[referencia]))
        lista_len.append(prueba)
    else:
        return lista_len
    

def matriz_uniq(referencia, dic):
    lista_len = []
    organismos = list(dic.keys())
    for org in organismos:
        comp = set(dic[org]) & set(dic[referencia])
        for org1 in organismos:
            if org1!=org and org1!=referencia:
                comp -= set(dic[org1])
        lista_len.append(len(comp))
    else:
        return lista_len
    

def uniq_lis(referencia, dic):
    dic_euk = {}
    organismos = list(dic.keys())
    for org in organismos:
        dic_euk[org] = []
        comp = set(dic[org]) & set(dic[referencia])
        for org1 in organismos:
            if org1!=org and org1!=referencia:
                comp -= set(dic[org1])
        dic_euk[org]=comp

    return dic_euk


def matr_cruzada(df,df1,title,file):
    lis = []
    lis1 = []
    for i in range(len(df)):
        lis.append(str(int(df.iloc[i,i])) + '\n' + str(int(df1.iloc[i,i])))
        lis1.append(int(df.iloc[i,i]))
        
    
    
    mat_0 = np.triu(df)
    mat_1 = np.tril(df1)
    mat_2 = mat_1 + mat_0
    df0 = pd.DataFrame(mat_2, columns=df.columns.tolist(), index=df.columns.tolist())
    df_t = df0.astype(int).astype(str)
    for i in range(len(df0)):
        df0.iloc[i,i] = int(lis1[i])
    for i in range(len(df_t)):
        df_t.iloc[i,i] = lis[i]
    t = np.array(df_t)
    
    fig, ax = plt.subplots(figsize=(25, 25))
    
    sb.heatmap(df0,robust=True,annot=t,fmt='',cmap='YlGnBu',annot_kws={'size':16},square=True,cbar_kws={"shrink": .7})
    ax.tick_params(labelsize=20)
    plt.rcParams['font.size'] = 20
    plt.title(title, loc='left',fontsize=20)
    plt.ylabel('')
    plt.savefig(file)
    plt.show()
    
def uniqueness(dic):
    clan_unique = {}
    tot = list(dic.keys())
    for group1 in tqdm(tot):
        tot = list(dic.keys())
        clan_unique[group1] = []
        tot.remove(group1)
        a = set(dic[group1])
        for group2 in tot:
            a -= set(dic[group2])
        clan_unique[group1] = a
    return clan_unique

def uniq_euk(dic):
    clan_unique = {}
    tot = list(dic.keys())
    tot.remove('Eukaryota')
    for group1 in tqdm(tot):
        tot = list(dic.keys())
        clan_unique[group1] = []
        tot.remove(group1)
        tot.remove('Eukaryota')
        a = set(dic[group1])
        for group2 in tot:
            a -= set(dic[group2])
        a = a & set(dic['Eukaryota'])
        clan_unique[group1] = a
    return clan_unique



#%%
############## FROM README UNIPROT REFERENCE PROTEOMES, REMOVE VIRUSES AND THEN DOWNLOAD EVERY PROTEOME
df=pd.DataFrame()
list_headers=[]
todo=[]
with open ("D:/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/README.txt",'r') as f:
    for linea in f:
        if linea.startswith('Proteome_ID'):
            headers=linea.split('\t')
            for x in headers:
                list_headers.append(x)
                #df.columns=list_headers
        if linea.startswith('UP0'):
            wlinea=linea.split('\t')
            todo.append(wlinea)
todo.insert(0,list_headers)


with open('D:/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/uniprot_dataset_junio.csv','w',newline='') as f2:
    writer = csv.writer(f2,delimiter=',')# quoting=csv.QUOTE_ALL,)  
    writer.writerow(list_headers) 
    writer.writerows(todo)     


df = pd.read_csv('D:/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/uniprot_dataset_junio.csv')

df = df[df.SUPERREGNUM!='viruses']

lis = df.Proteome_ID.tolist()

# lis = ['UP000289406','UP000320513']
urls = []
for i in lis:
    urls.append((i+'.fasta',"https://www.uniprot.org/uniprot/?query=proteome:" + i + "&format=fasta"))


results = ThreadPool(8).imap_unordered(fetch_url, urls)
for path in tqdm(results):
    print(path)



#%%
"""
########################## FROM PFAM_SCAN TO 1 raw table with every architecture
"""

path='D:/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/'
database='D:/SINCRONIZACION_CABD/new_tax_red/'
protein_comb_new_NEWDATA_repeats(path,'hmm_acc','clan',database,'definitive09112022')

#%%

############################PLOTS
  
infile = '/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/todas_arq_comb_definitive09112022.csv'

# infile = '/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/todas_arq_comb_repeats_info_FULL.csv'

with open (infile,'r') as f:
    org_king = {}
    org_phyl = {}
    org_prot = {}
    dic_multi_arch = {}
    dic_uni_arch = {}
    dic_bi = {}
    dic_unit = {}
    dic_comb = {}
    dic_rep = {}
    for line in tqdm(f):
        wline = line.rstrip('\n').split(',')
        arch = wline[6].split('-')
        if wline[2] not in org_king:
            org_king[wline[2]] = []
            org_phyl[wline[2]] = []
            org_prot[wline[2]] = []
            dic_multi_arch[wline[2]] = []
            dic_uni_arch[wline[2]] = []
            dic_bi[wline[2]] = []
            dic_unit[wline[2]] = []
            dic_comb[wline[2]] = []
            dic_rep[wline[2]] = []
        if wline[4] not in org_king[wline[2]]:
            org_king[wline[2]]=wline[4]
        if wline[3] not in org_phyl[wline[2]]:
            org_phyl[wline[2]]=wline[3]
        if wline[0] not in org_phyl[wline[2]] and len(wline[6].split('-'))>1: #CHEQUEAR!
            org_prot[wline[2]].append(wline[0])
        dic_comb[wline[2]].append(wline[6])
        if 'Repeat' in wline[8]:
            dic_rep[wline[2]].append(wline[6])
        for i in arch:
        # for i,x in zip(arch, wline[8].split('-')):
            # if x=='Repeat':
                dic_unit[wline[2]].append(i)
        for i in range(0,len(arch)-1):
            # if 'Repeat' in wline[8]:
            bigram = arch[i] + '-' + arch[i+1]
            dic_bi[wline[2]].append(bigram)
        if len(wline[6].split('-'))>1 and not 'Repeat' in wline[8]:
            dic_multi_arch[wline[2]].append(wline[6])
        if len(wline[6].split('-'))==1 and not 'Repeat' in wline[8]:
            dic_uni_arch[wline[2]].append(wline[6])


for i in dic_uni_arch:
    dic_uni_arch[i] = list(set(dic_uni_arch[i]))
    dic_multi_arch[i] = list(set(dic_multi_arch[i]))    
    dic_bi[i] = list(set(dic_bi[i]))
    dic_unit[i] = list(set(dic_unit[i]))
    dic_comb[i] = list(set(dic_comb[i]))
    dic_rep[i] = list(set(dic_rep[i]))
   


lis_tot = []
lis_tot.append(dic_multi_arch)
lis_tot.append(dic_uni_arch)
lis_tot.append(dic_bi)
lis_tot.append(dic_unit)
lis_tot.append(dic_comb)
lis_tot.append(dic_rep)
lis_tot.append(org_prot)
lis_tot1 = ['dic_multi_arch', 'dic_uni_arch', 'dic_bi', 'dic_unit','dic_comb','dic_rep','org_prot']



def get_phyla_TOT(dic):
    dic_phyl = {}
    for i in dic:
        phyl = org_phyl[i]
        if phyl not in dic_phyl:
            dic_phyl[phyl] = []
        for x in dic[i]:
            dic_phyl[phyl].append(x)
    for i in dic_phyl:   
        dic_phyl[i] = list(set(dic_phyl[i]))
    return dic_phyl


lis_tot_phyl = []
lis_tot_phyl1 = []
for i,j in zip(lis_tot,lis_tot1):
    exec(f'{j}_phyl = get_phyla_TOT(i)')
    lis_tot_phyl1.append(f'{j}_phyl')
    lis_tot_phyl.append(get_phyla_TOT(i))



dic_bi_uniq = uniqueness(dic_bi_phyl)
dic_comb_uniq = uniqueness(dic_comb_phyl)
dic_unit_uniq = uniqueness(dic_unit_phyl)
dic_uni_arch_uniq = uniqueness(dic_uni_arch_phyl)
dic_multi_arch_uniq = uniqueness(dic_multi_arch_phyl)
dic_rep_uniq = uniqueness(dic_rep_phyl)
    
import math
df = pd.DataFrame(index=list(dic_unit.keys()))

# df['org'] = ''
df['superkingdom'] = int()
df['phylum'] = int()
# df['unigrams_typ'] = int()
df['unigrams'] = int()
# df['bigrams_typ'] = int()
# df['unigrams_uniqo'] = int()
df['bigrams'] = int()
# df['bigrams_uniqo'] = int()
# df['multi_arch_typ'] = int()
df['multi_arch'] = int()
# df['multi_arch_uniqo'] = int()
# df['uni_arch_typ'] = int()
df['uni_arch'] = int()
# df['uni_arch_uniqo'] = int()
# df['total_arch_typ'] = int()
df['total_arch'] = int()
# df['total_arch_uniqo'] = int()
# df['rep_arch_typ'] = int()
df['rep_arch'] = int()
# df['rep_arch_uniqo'] = int()
df['complex'] = int()
# df['prop_multi'] = int()
# df['n_protein'] = int()
for i in df.index:
    # df.org[i] = i
    df['superkingdom'][i] = org_king[i]
    df['phylum'][i] = org_phyl[i]
    # df['unigrams_typ'][i] = int(len(dic_unit_typ[i]))
    df['unigrams'][i] = int(len(dic_unit[i]))
    # df['unigrams_uniqo'][i] = int(len(dic_unit_uniqo[i]))
    df['bigrams'][i] = int(len(dic_bi[i]))
    # df['bigrams_typ'][i] = int(len(dic_bi_typ[i]))
    # df['bigrams_uniqo'][i] = int(len(dic_bi_uniqo[i]))
    # df['multi_arch_typ'][i] = int(len(dic_multi_arch_typ[i]))
    df['multi_arch'][i] = int(len(dic_multi_arch[i]))
    # df['multi_arch_uniqo'][i] = int(len(dic_multi_arch_uniqo[i]))
    # df['uni_arch_typ'][i] = int(len(dic_uni_arch_typ[i]))
    df['uni_arch'][i] = int(len(dic_uni_arch[i]))
    # df['uni_arch_uniqo'][i] = int(len(dic_uni_arch_uniqo[i]))
    # df['total_arch'][i] =  int(len(dic_comb[i]))
    # df['rep_arch'][i] =  int(len(dic_rep[i]))
    # df['total_arch_typ'][i] =  len(dic_comb_typ[i])
    df['total_arch'][i] =  len(dic_comb[i])
    # df['total_arch_uniqo'][i] =  len(dic_comb_uniqo[i])
    # df['rep_arch_typ'][i] =  len(dic_rep_typ[i])
    df['rep_arch'][i] =  len(dic_rep[i])
    # df['rep_arch_uniqo'][i] =  len(dic_rep_uniqo[i])
    df['complex'][i] = len(dic_rep[i])+len(dic_multi_arch[i])
    # df['prop_multi'][i] = len(dic_multi_arch[i])/len(dic_comb[i])
    # df['n_protein'][i] = len(org_prot[i])


dic_count = df.phylum.value_counts()
dic_count[i]
# for i in df1.index:
#     dic_count[i] = df.phylum.value_counts()

for i in df.index:
    print(i)
    # print(i, str(df['phylum'].value_counts(i)))
    df['phylum'][i] = df['phylum'][i] + ' (n=' + str(dic_count[df.phylum[i]]) + ')'
                                   
# df.phylum.value_counts()           
                                   
# df['phylum'] = df[]
# str(df['phylum'].value_counts()[df.phylum[df.index==i]])

df1 = pd.DataFrame(index=list(dic_rep_phyl.keys()))
df1['phylum'] = str()
df1['unigrams'] = int()
df1['unigrams_uniq'] = int()
df1['bigrams'] = int()
df1['bigrams_uniq'] = int()
df1['total_arch'] = int()
df1['total_arch_uniq'] = int()
df1['multi_arch'] = int()
df1['multi_arch_uniq'] = int()
df1['uni_arch'] = int()
df1['uni_arch_uniq'] = int()
df1['rep_arch'] = int()
df1['rep_arch_uniq'] = int()
for i in df1.index:
    # df1['unigrams'][i] = (len(dic_unit_phyl[i])-len(dic_unit_uniq[i]))
    df1['phylum'][i] = i + ' (n=' + str(dic_count[i]) + ')'
    df1['unigrams'][i] = len(dic_unit_phyl[i])
    df1['unigrams_uniq'][i] = len(dic_unit_uniq[i])
    # df1['bigrams'][i] = (len(dic_bi_phyl[i])-len(dic_bi_uniq[i]))
    df1['bigrams'][i] = len(dic_bi_phyl[i])
    df1['bigrams_uniq'][i] = len(dic_bi_uniq[i])
    df1['total_arch'][i] = len(dic_comb_phyl[i])
    df1['total_arch_uniq'][i] = len(dic_comb_uniq[i])
    df1['multi_arch'][i] = len(dic_multi_arch_phyl[i])
    df1['multi_arch_uniq'][i] = len(dic_multi_arch_uniq[i])
    df1['uni_arch'][i] = len(dic_uni_arch_phyl[i])
    df1['uni_arch_uniq'][i] = len(dic_uni_arch_uniq[i])
    df1['rep_arch'][i] = len(dic_rep_phyl[i])
    df1['rep_arch_uniq'][i] = len(dic_rep_uniq[i])
    

import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))

fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['unigrams'],
    x=df.phylum,# + ' (n=' + str(df['phylum'].value_counts()[df1.index]) + ')',
    name='Uni-grams',
    marker_color='cadetblue',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=df1['unigrams']-df1['unigrams_uniq'],
    x=df1.phylum,# + ' (n=' + str(df['phylum'].value_counts()[df1.index]) + ')',
    name='Cum. Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.add_trace(go.Bar(
    y=df1['unigrams_uniq'],
    x=df1.phylum,# + ' (n=' + str(df['phylum'].value_counts()[df1.index]) + ')',
    name='Cum. Uni-grams Uniq',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.update_layout(barmode='stack',title="Uni-grams",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_unigrams.html")

# df['phylum'].value_counts()['Halobacteriota']


import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))

fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['bigrams']/df['unigrams'],
    x=df.phylum,
    name='Bi-grams/Uni-grams',
    marker_color='chocolate',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=(df1['bigrams']-df1['bigrams_uniq'])/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Bi-grams/Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.add_trace(go.Bar(
    y=df1['bigrams_uniq']/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Bi-grams Uniq/Uni-grams',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)#,secondary_y=False)
# fig.update_layout(barmode='stack',title="Bigrams")
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_bigrams-unigrams.html")
fig.update_layout(barmode='stack',title="Bi-grams",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_bigrams.html")





import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))


fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['total_arch']/df['unigrams'],
    x=df.phylum,
    name='Total Arch/Uni-grams',
    marker_color='violet',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=(df1['total_arch']-df1['total_arch_uniq'])/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Total Arch/Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.add_trace(go.Bar(
    y=df1['total_arch_uniq']/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Total Arch Uniq/Uni-grams',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)#,secondary_y=False)
# fig.update_layout(barmode='stack',title="Total Architectures")
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_totalarch-unigrams.html")
fig.update_layout(barmode='stack',title="Total Architectures",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_total_arch.html")


import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))


fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['uni_arch']/df['unigrams'],
    x=df.phylum,
    name='Single-Domain Arch/Uni-grams',
    marker_color='yellowgreen',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=(df1['uni_arch']-df1['uni_arch_uniq'])/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Single-Domain Arch/Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.add_trace(go.Bar(
    y=df1['uni_arch_uniq']/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Single-Domain Arch Uniq/Uni-grams',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)#,secondary_y=False)
# fig.update_layout(barmode='stack',title="Single-Domain Architectures")
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_singlearch-unigrams.html")
fig.update_layout(barmode='stack',title="Single-Domain Architectures",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_single_dom_arch.html")


import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))


fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['multi_arch']/df['unigrams'],
    x=df.phylum,
    name='Multi Arch/Uni-grams',
    marker_color='cornflowerblue',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=(df1['multi_arch']-df1['multi_arch_uniq'])/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Multi Arch/Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)#,secondary_y=False)
fig.add_trace(go.Bar(
    y=df1['multi_arch_uniq']/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Multi Arch Uniq/Uni-grams',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)#,secondary_y=False)
# fig.update_layout(barmode='stack',title="Multi-Domain Architectures")
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_multiarch-unigrams.html")
fig.update_layout(barmode='stack',title="Multi-Domain Architectures",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_Multi-Domain_arch.html")


import plotly.graph_objects as go
from plotly.subplots import make_subplots
fig = make_subplots(shared_xaxes=True,rows=2, cols=1,specs=[[{"secondary_y": True},],[{"secondary_y": True}]],subplot_titles=("A)","B)"))


fig.add_trace(go.Box(
    # defining y axis in corresponding
    # to x-axis
    y=df['rep_arch']/df['unigrams'],
    x=df.phylum,
    name='Rep Arch/Uni-grams',
    marker_color='burlywood',hovertext=df.index, boxpoints='all'),row=1,col=1)#,secondary_y=False)#,name="Proportion of multi-domain architectures")
fig.add_trace(go.Bar(
    y=(df1['rep_arch']-df1['rep_arch_uniq'])/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Rep Arch/Uni-grams',
    marker_color='gold',opacity=0.4),row=2,col=1)
fig.add_trace(go.Bar(
    y=df1['rep_arch_uniq']/df1['unigrams'],
    x=df1.phylum,
    name='Cum. Rep Arch Uniq/Uni-grams',
    marker_color='blueviolet',opacity=0.4),row=2,col=1)

# fig.update_layout(barmode='stack',title="Repeat Architectures")#,
                  # title = 'Stacked bar chart!',
                  # template = 'plotly_dark').show()
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_reparch-unigrams.html")
fig.update_layout(barmode='stack',title="Repeat Architectures",font=dict(size=30))
fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_xaxes(tickangle=50)

for i in fig['layout']['annotations']:
    i['font'] = dict(size=25)
        # family="Courier New, monospace",
        
        # color="RebeccaPurple"))
# plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/" + "box_unigrams.html")
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "box_Repeat_arch.html")

#####BAR PLOT PROPORTION

df2 = pd.DataFrame(index=list(dic_rep_phyl.keys()))
df2['phylum'] = str()
df2['uni_arch'] = int()
df2['rep_arch'] = int()
df2['multi_arch'] = int()

df2['complex'] = int()

for i in df2.index:
    df2['phylum'][i] = i + ' (n=' + str(dic_count[i]) + ')'
    df2['rep_arch'][i] = len(dic_rep_phyl[i])/len(dic_comb_phyl[i])
    df2['uni_arch'][i] = len(dic_uni_arch_phyl[i])/len(dic_comb_phyl[i])
    df2['multi_arch'][i] = len(dic_multi_arch_phyl[i])/len(dic_comb_phyl[i])
    df2['complex'][i] = (len(dic_multi_arch_phyl[i])+len(dic_rep_phyl[i]))/len(dic_comb_phyl[i])
    

df2.sort_values(by = ['complex'],ascending = True ,inplace=True)  

ax = df2.iloc[:,0:4].plot(kind='barh', 
         stacked=True, 
         colormap='Pastel2', 
         figsize=(8, 8),alpha=0.8,width=0.9)
plt.legend(
        bbox_to_anchor=(0.5, 1.02),
        loc="lower center",
        borderaxespad=0,
        frameon=False,
        ncol=3,
    )

plt.xlabel("Proportion of architectures")
plt.ylabel("Phyla")
plt.tight_layout()
ax.set_yticklabels(df2.phylum)
plt.show()
