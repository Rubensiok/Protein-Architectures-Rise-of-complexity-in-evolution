#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 19:16:38 2023

@author: ruben
"""

import pandas as pd
from tqdm import tqdm
import plotly
from matplotlib import pyplot as plt
# %matplotlib inline
import plotly.express as px
from goatools.godag.prtfncs import GoeaPrintFunctions
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
import goatools

a = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/libreraria_GO_all_protein_pruebass.csv',header=0,index_col=0)
assoc = {}
for i in tqdm(a.index):
    if i not in assoc:
        assoc[i] = set()
    assoc[i].add(str(a['GO_ID'][i]))
del a
go = obo_parser.GODag('/home/ruben/SINCRONIZACION_CABD/bigrams/go-basic.obo')

grac = ['Proteobacteria','Spirochaetota','Bacteroidota','Acidobacteriota',
            'Planctomycetota','Verrucomicrobiota','Chlamydiota','Bdellovibrionota','Campylobacterota','Desulfobacterota']
terra = ['Actinobacteriota','Firmicutes','Cyanobacteria','Chloroflexota',
         'Aquificota','Deinococcota','Synergistota','Fusobacteriota','Thermotogota']

#%%

# 1 per round

df = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_tfidf_EUK_PROK1.csv',index_col = 0)

df =pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_TDIDF_Archaea_Bacteria_1.csv', index_col = 0)

df = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/uniarch_and_bigrams_tfidf_grac_terra_1.csv',index_col = 0)

df = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_tfidf_Plancto_GRAC1.csv',index_col = 0)

df = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_TDIDF_Cyano_vs_Terra_1.csv',index_col = 0)

df = pd.read_csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_TDIDF_PROTEO_VS_GRACI1.csv',index_col = 0)




#%%

df.index = df.index.str.replace('.','-')
dic_contr = {}
for i in df.GroupContrib:
    dic_contr[i] = []
    
for i in df.index:
    dic_contr[df.GroupContrib[df.index==i].iloc[0]].append(i.replace('.','-'))
    
infile = '/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/todas_arq_comb_definitive09112022.csv'

dic_acc = {}
for i in dic_contr:
    for x in dic_contr[i]:
       dic_acc[x] = []
    # dic_acc[dic_contr[]] = []

#%%
  
# Select depends on the target

with open(infile) as f:
    for line in tqdm(f):
        wline = line.split(',')
        for i in dic_contr:
            if i==wline[4]:
                for x in dic_contr[i]:
                    if len(x.split('-'))==2 and x in wline[6]:
                        dic_acc[x].append(wline[0])
                        # print(wline[0],wline[8])
                    if len(x.split('-'))==1 and x==wline[6]:
                        dic_acc[x].append(wline[0])
                        # print(wline[0],wline[8])


with open(infile) as f:
    for line in tqdm(f):
        wline = line.split(',')
        # if wline[3] in grac:
        for i in dic_contr:
            for x in dic_contr[i]:
                if (i=='Gracillicutes' and wline[3] in grac) or (i=='Terrabacteria' and wline[3] in terra):
                    if len(x.split('-'))==2 and x in wline[6]:
                        dic_acc[x].append(wline[0])
                    if len(x.split('-'))==1 and x==wline[6]:
                        dic_acc[x].append(wline[0])

with open(infile) as f:
    for line in tqdm(f):
        wline = line.split(',')
        for i in dic_contr:
            if i==wline[3]:
                for x in dic_contr[i]:
                    if len(x.split('-'))==2 and x in wline[6]:
                        dic_acc[x].append(wline[0])
                        # print(wline[0],wline[8])
                    if len(x.split('-'))==1 and x==wline[6]:
                        dic_acc[x].append(wline[0])
                        # print(wline[0],wline[8])

#%%

for i in dic_acc:
    dic_acc[i] = list(set(dic_acc[i]))
    
dic_go = {}
for i in dic_acc:
    dic_go[i] = []
lis_go = []
for i in dic_go:
    for x in dic_acc[i]:
        try:
            dic_go[i].append(','.join(assoc[x]))
            lis_go.append(','.join(assoc[x]))
        except:
            continue
lis_go = list(set(lis_go))

df1 = pd.DataFrame(index=dic_go.keys())
df1['goes'] = ''
for i in df1.index:
    df1['goes'][i] = ';'.join(list(set(dic_go[i])))
df1['lis_go'] = ''
for i in df1.goes:
    # print(i)
    i1 = i.split(';')
    print(i1)
    c = []
    for x in i1:
        if x!='':
            try:
                b = x + ' ' + go[x].name
            except:
                b = x + ' NO_GO_name'
                # pass
        else:
            b = 'NO_GO_term'
        c.append(b)
    df1.lis_go[df1.goes==i] = ';'.join(c)

df1['count_go'] = ''
for i in df1.index:
    i1 = df1.lis_go[i].split(';')
    c = []
    for x in i1:
        x1 = x.split(' ')
        b = dic_go[i].count(x1[0])
        c.append(str(b))
    df1.count_go[i] = ';'.join(c)
        
dic_count = {}
dic_fin = {}
for i in dic_go:
    dic_count[i] = int()
    dic_fin[i] = []
    for x in lis_go:
        a = dic_go[i].count(x)
        # dic_count[i].append(str(a) + ' ' + x)
        if a>dic_count[i]:
            dic_count[i] = a
            dic_fin[i] = x
            print(i,x,a)

df1['splsda_importance'] = int()
for i in df1.index:
    df1['splsda_importance'][i] = df.importance[df.index==i].iloc[0]

df1['GroupContrib'] = ''
for i in df1.index:
    df1['GroupContrib'][i] = df.GroupContrib[df.index==i].iloc[0]

df1['plotted'] = ''
for i,j in zip(df1.lis_go,df1.count_go):
    i1 = i.split(';')
    j1 = j.split(';')
    j1 = [int(x) for x in j1]
    # for x in j1:
    #     print(type(x))
    c = [z for _,z in sorted(zip(j1,i1),reverse=True)]
    c = c[0:3]
    print(c)
    df1.plotted[df1.lis_go==i] = ';'.join(c)
    
#%%

# Select depend on the target

import plotly.express as px
import plotly

df1.splsda_importance = abs(df1.splsda_importance)
df3 = df1[df1.GroupContrib=='Eukaryota']
df3 = df3.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Eukaryota functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=18,textposition="auto",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_EUK_PROKA" + ".html")


df1.splsda_importance = abs(df1.splsda_importance)
# df3 = df1[df1.GroupContrib=='Archaea']
df3 = df1.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Archaea functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=18,textposition="inside",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_arch_bac" + ".html")



df1.splsda_importance = abs(df1.splsda_importance)
# df3 = df1[df1.GroupContrib=='Gracillicutes']
df3 = df1.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Gracillicutes/Terrabacteria functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=17,textposition="auto",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_gracill_terra" + ".html")



df1.splsda_importance = abs(df1.splsda_importance)
# df3 = df1[df1.GroupContrib=='Gracillicutes']
df3 = df1.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Planctomycetota functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=18,textposition="auto",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_plancto_grac" + ".html")



df1.splsda_importance = abs(df1.splsda_importance)
# df3 = df1[df1.GroupContrib=='Gracillicutes']
df3 = df1.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Cyanobacteria functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=16,textposition="auto",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_Cyano_terra" + ".html")



df1.splsda_importance = abs(df1.splsda_importance)
# df3 = df1[df1.GroupContrib=='Gracillicutes']
df3 = df1.sort_values(by='splsda_importance', ascending=False)
df3 = df3.iloc[0:15,:]
df3 = df3.sort_values(by='splsda_importance', ascending=    True)
fig = px.bar(df3, x=df3.splsda_importance, y=df3.index, text=df3.plotted,color = df3.GroupContrib,orientation = 'h',title = 'Proteobacteria functional Analysis',labels={
                     "index": "Bi-grams/Single Domain Architectures"})
                
fig.show()  
fig.update_traces(textfont_size=18,textposition="auto",cliponaxis=True,insidetextanchor = 'start')
plotly.offline.plot(fig, filename= "/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/final_plots/" + "bar_bi_arch_proteo_grac" + ".html")