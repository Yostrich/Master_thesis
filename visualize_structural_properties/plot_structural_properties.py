# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:43:50 2023

@author: Danish
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import ks_2samp

ident="pooled"

# Read the CSV files


# PAK_P3
#TSS_class_PAKP3 = [0]*3 +[1]*2 + [0]*2 +[1]*1 + [0]*4 + [1]*3 + [0]*1 +[1]*2 + [0]*1 +[1]*1 + [0]*2 +[1]*2 + [0]*1 +[1]*2 + [0]*1 +[1]*5 + [0]*2 +[1]*1+ [0]*1  +[1]*1 + [0]*1 +[1]*1 +[0]*1 +[1]*1 + [0]*1 +[1]*3 + [0]*1 +[1]*4 + [0]*1 +[1]*3 + [0]*2 +[1]*3 + [0]*1 +[1]*1 + [0]*2 +[1]*1 + [0]*1 + [1]*3 + [0]*1 + [1]*4
#TTS_class_PAKP3 = [0]*5 +[1]*5 + [0]*5 +[1]*2 + [0]*2 +[1]*1 + [0]*4 +[1]*1 + [0]*1 +[1]*1 + [0]*1 +[1]*1 +[0]*1
TSS_class_PAKP3 = [0]*74
TTS_class_PAKP3 = [0]*30

#LUZ19 motifs
#TSS_class_LUZ19 = [1]*2 +[0]*1 + [1]*2 + [0]*4 +[2]*3
#TTS_class_LUZ19 = [0]*7 + [1] + [0] + [0]*3

TSS_class_LUZ19 = [0]*12
TTS_class_LUZ19 = [0]*12

#LUZ100
#TSS_class_LUZ100 = [1]*3 + [2]*3 + [0]*1 + [2]*1 + [0]*2 + [2]*2
#TTS_class_LUZ100 = [0]*2 + [1] + [0] + [1]*2 + [0] + [0] + [0]*4 + [1] + [0]
TSS_class_LUZ100 = [0]*12
TTS_class_LUZ100 = [0]*14


#14-1
#TSS_class_fto = [1]*15 + [0]*2 + [1]
#TTS_class_fto = [0] + [0] + [0]*2 + [1]*4 + [0]*2 + [1] + [0]
TSS_class_fto = [0]*18
TTS_class_fto = [0]*12



#phiKZ
#TSS_class_phiKZ = [1]*2 + [0]*2 + [2]*1 + [0]*5 + [1]*1 + [3]*1 + [0]*2 +[3]*1 + [0]*1 + [1]*1 + [2]*1 + [3]*1 + [0]*1 + [3]*1 + [0]*1 + [3]*1 + [2]*1 + [1]*2 + [3]*1 + [1]*1 + [0]*1 + [2]*1 + [1]*1 + [0]*2 + [2]*1 + [3]*1+ [0]*1 + [3]*1 + [1]*1 + [0]*2 + [1]*1 + [0]*2 + [3]*2 + [0]*2 + [2]*1 + [0]*3 + [2]*1 + [0]*1 +[2]*1 + [0]*1 + [3]*1 + [0]*2 + [2]*3 + [0]*1 + [2]*1 + [1]*2 + [1]*2 + [0]*5 + [3]*1 + [0]*5 + [2]*1 + [1]*1 +[0]*1 + [2]*1 + [0]*10 + [1]*2 + [0]*1 + [3]*1 + [2]*1 + [0]*2 + [1]*1 + [0]*1 + [2]*1 + [0]*3 + [3]*1 + [0]*1 + [1]*1 + [0]*1 + [3]*1 + [1]*1 + [3]*1 + [0]*2 + [2]*2 + [1]*1 + [0]*1 + [3]*1 + [0]*1 + [1]*2 + [3]*1 + [1]*1 + [0]*1 + [3]*1 + [1]*1 + [3]*1 + [0]*4 + [2]*1 + [3]*1 + [2]*1 + [1]*1 + [0] + [0]*4 + [2]*3 + [3]*1 + [2]*1 + [0]*2  + [1]*2 + [0]*3 +[3]*1 +[2]*1 + [0]*6 + [3]*1 + [0]*4 + [3]*1 + [0]*1 + [2]*1 + [3]*1 + [0]*1 + [3]*1 + [0]*2 + [2]*1 + [0]*5 + [2]*1 + [3]*1 + [2]*1 + [0]*4 + [3]*1 + [0]*1 + [2]*3 + [3]*7 + [0]*1 + [3]*2 + [0]*2 + [2]*1 + [3]*1 + [0]*1
#TTS_class_phiKZ = [1]*1 + [0]*2 + [1]*2 + [0]*3 + [1]*1 + [0]*1 + [1]*1 + [0]*2 + [1]*1 + [0]*8 + [1]*1 + [0]*2 + [1]*2 + [0]*3 + [1]*1 + [0]*1 + [1]*1 + [0]*1+ [1]*1 + [0]*1+ [1]*1 + [0]*1+ [1]*1 + [0]*1+ [1]*1 + [0]*3 + [1]*2 + [0]*2 + [1]*2 + [0]*2 + [1]*2 + [0]*1 + [1]*10 + [0]*3 + [1]*4 + [0]*1 + [1]*1 + [0]*1 + [1]*2 + [0]*1 + [1]*2 + [0]*5 + [1]*2 + [0]*1 + [1]*4 + [0]*1 + [1]*4 + [0]*1 + [1]*1 + [0]*1 + [1]*3
TSS_class_phiKZ = [0]*207
TTS_class_phiKZ = [0]*103


#LUZ7
#TSS_class_LUZ7 = [0]*5 + [1]*5 + [2]*76
#TTS_class_LUZ7 =  [0]*3 + [1] + [0]*8 + [0] + [0]*5 + [1] + [0]*3 + [1] + [0]*4 + [1] + [0]*2 + [1] + [0]*8 + [1] + [0]*9 + [1] + [0]*2 + [1] 
TSS_class_LUZ7 = [0]*86
TTS_class_LUZ7 = [0]*53

list_phages = {'LUZ7':[TSS_class_LUZ7,TTS_class_LUZ7], 'LUZ100':[TSS_class_LUZ100,TTS_class_LUZ100],
               'LUZ19':[TSS_class_LUZ19,TTS_class_LUZ19], 'phiKZ':[TSS_class_phiKZ,TTS_class_phiKZ],
               '14-1':[TSS_class_fto,TTS_class_fto],'PAKP3':[TSS_class_PAKP3,TTS_class_PAKP3]}


def plotAverageStructuralPropAndDens(y_values_RE, lower, upper, classList, colorList, y_values_rand, titleName, yName, xName, TXS, filename):
    x_axis = y_values_RE[lower:upper].iloc[:,0]
    y_values_RE = y_values_RE[lower:upper].iloc[:,1:].T
    y_values_rand = y_values_rand[lower:upper].iloc[:,1:]
    y_values_RE['category'] = classList
    groups = y_values_RE.groupby('category')
    
    #data = pd.DataFrame({'Random': y_values_rand.mean(axis=1)})
    data = {'Random': y_values_rand.mean(axis=0)}
    
    for name, group in groups:
        test = group.iloc[:,0:-1].T
        plt.plot(x_axis, test.mean(axis = 1), color = colorList[name])
        #axis[1].fill_between(x_axis,test.mean(axis=1) - test.std(axis=1), test.mean(axis=1) + test.std(axis=1), color = colorList[name], alpha=0.2)
        data[name] = test.mean(axis = 0)
        ks_test = ks_2samp(data['Random'], test.mean(axis = 0))
        print('ks-test {prop} (group: {group}): pval = {pval}'.format(prop = titleName, pval = ks_test[1], group = name))
    
    #axis[0].boxplot(data, labels = boxplotNames)
    plt.plot(x_axis, y_values_rand.mean(axis=1), color = "black")
    plt.title("{name}".format(name = titleName))
    plt.ylabel(yName)
    plt.xlabel(xName)
    xaxis_vis = list(filter(lambda x: x % 10 == 0 ,x_axis))
    ticks_vis = [ i if i != 0 else TXS for i in xaxis_vis]
    plt.xticks(xaxis_vis,ticks_vis)
    plt.savefig('{fileName}_{T}_profile.png'.format(fileName = filename, T = TXS))
    plt.close()
    keys = list(data.keys())
    vals = [data[k][:-1] for k in keys]
    for i in range(0,len(vals)):
        sns.kdeplot(vals[i], color = (["black"] + colorList)[i])
    plt.xlabel(yName)
    plt.title("{name}".format(name = titleName))
    plt.savefig('results/{fileName}_{T}_dist.png'.format(fileName = filename, T = TXS))
    plt.close()
    
def plotAverageStructralPropertiesAllPhages(listPhages, DNAprop, lower, upper, TXS, y_axis):
    for phage in listPhages:
        yval = pd.read_csv('structural_properties_data/{T}_{phage}_pooled_{prop}.csv'.format(phage=phage, prop =  DNAprop, T = TXS))
        x_axis = yval[lower:upper].iloc[:,0]
        y_values_RE = yval[lower:upper].iloc[:,1:]
        plt.plot(x_axis, y_values_RE.mean(axis=1))
    plt.legend(listPhages)
    plt.xlabel("position relative to TSS")
    plt.ylabel(y_axis)
    xaxis_vis = list(filter(lambda x: x % 10 == 0 ,x_axis))
    ticks_vis = [ i if i != 0 else TXS for i in xaxis_vis]
    plt.xticks(xaxis_vis,ticks_vis)
    plt.savefig('results/all_phages_{T}_{prop}.png'.format(prop =  DNAprop, T = TXS))
    plt.close()
    
for phage in  list_phages:
    TSS_curv = pd.read_csv('structural_properties_data/TSS_{phage}_{ident}_curvature.csv'.format(phage=phage, ident=ident))
    TSS_bend = pd.read_csv('structural_properties_data/TSS_{phage}_{ident}_bendability.csv'.format(phage=phage, ident=ident))
    TSS_stab = pd.read_csv('structural_properties_data/TSS_{phage}_{ident}_stability.csv'.format(phage=phage, ident=ident))

    TTS_curv = pd.read_csv('structural_properties_data/TTS_{phage}_{ident}_curvature.csv'.format(phage=phage, ident=ident))
    TTS_bend = pd.read_csv('structural_properties_data/TTS_{phage}_{ident}_bendability.csv'.format(phage=phage, ident=ident))
    TTS_stab = pd.read_csv('structural_properties_data/TTS_{phage}_{ident}_stability.csv'.format(phage=phage, ident=ident))

    TSS_curv_rand = pd.read_csv('structural_properties_data/TSS_random_seq_{phage}_{ident}_curvature.csv'.format(phage=phage, ident=ident))
    TSS_bend_rand = pd.read_csv('structural_properties_data/TSS_random_seq_{phage}_{ident}_bendability.csv'.format(phage=phage, ident=ident))
    TSS_stab_rand = pd.read_csv('structural_properties_data/TSS_random_seq_{phage}_{ident}_stability.csv'.format(phage=phage, ident=ident))

    TTS_curv_rand = pd.read_csv('structural_properties_data/TTS_random_seq_{phage}_{ident}_curvature.csv'.format(phage=phage, ident=ident))
    TTS_bend_rand = pd.read_csv('structural_properties_data/TTS_random_seq_{phage}_{ident}_bendability.csv'.format(phage=phage, ident=ident))
    TTS_stab_rand = pd.read_csv('structural_properties_data/TTS_random_seq_{phage}_{ident}_stability.csv'.format(phage=phage, ident=ident))
    
    TSS_class = list_phages[phage][0]
    TTS_class = list_phages[phage][1]
    plotAverageStructuralPropAndDens(TSS_curv, 20, 81 ,TSS_class, ["red","blue", "green", 'orange'], TSS_curv_rand, "{phage}".format(phage=phage), "Curvature (d/Imax)", "position relative to TSS", 'TSS', "{phage}_curvature".format(phage=phage))
    plotAverageStructuralPropAndDens(TSS_stab, 20, 81 ,TSS_class, ["red","blue", "green", 'orange'], TSS_stab_rand, "{phage}".format(phage=phage), "AFE (kcal/mol)", "position relative to TSS", 'TSS', "{phage}_stability".format(phage=phage))
    plotAverageStructuralPropAndDens(TSS_bend, 20, 81 ,TSS_class, ["red","blue", "green", 'orange'], TSS_bend_rand, "{phage}".format(phage=phage), "Bendability", "position relative to TSS", 'TSS', "{phage}_bendability".format(phage=phage))
    plotAverageStructuralPropAndDens(TTS_curv, 20, 81 ,TTS_class, ["red","blue", "green", 'orange'], TSS_curv_rand, "{phage}".format(phage=phage), "Curvature (d/Imax)", "position relative to TTS", 'TTS', "{phage}_curvature".format(phage=phage))
    plotAverageStructuralPropAndDens(TTS_stab, 20, 81 ,TTS_class, ["red","blue", "green", 'orange'], TSS_stab_rand, "{phage}".format(phage=phage), "AFE (kcal/mol)", "position relative to TTS", 'TTS', "{phage}_stability".format(phage=phage))
    plotAverageStructuralPropAndDens(TTS_bend, 20, 81 ,TTS_class, ["red","blue", "green", 'orange'], TSS_bend_rand,"{phage}".format(phage=phage), "Bendability", "position relative to TTS", 'TTS', "{phage}_bendability".format(phage=phage))


plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'curvature', 20, 81, 'TSS', "Curvature (d/Imax)")
plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'bendability', 20, 81, 'TSS', 'Bendability')
plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'stability', 20, 81, 'TSS', "AFE (kcal/mol)")
plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'curvature', 20, 81, 'TTS', "Curvature (d/Imax)")
plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'bendability', 20, 81, 'TTS', 'Bendability')
plotAverageStructralPropertiesAllPhages(['LUZ7', 'LUZ100', 'LUZ19', 'phiKZ', '14-1','PAKP3'], 'stability', 20, 81, 'TTS', "AFE (kcal/mol)")
