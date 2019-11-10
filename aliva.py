#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : AliVa
# File Name       : aliva.py
# Version       : 1.0.3
# Encoding        : python
# Creation Date   : 2019/09/1
# Author : Taro Maeda 
# license     MIT License (http://opensource.org/licenses/mit-license.php)
# Copyright (c) 2019 Taro Maeda
# ======================================================================
# 



import argparse
import pandas as pd
from Bio import SeqIO
from Bio import Seq

import matplotlib.pyplot as plt
import seaborn as sns; sns.set() #seabornライブラリを読み込み、スタイルをセットする
import numpy as np
import pandas as pd



def GET_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="File path to a genome sequence file", required=True)
    parser.add_argument('-w','--win', help="window size (bp) (default = 50)", required=False)
    parser.add_argument('-p','--plot', help="Manualy modify the plot with seaborn window", action='store_true', required=False)
    parser.add_argument('-o','--out', help="Prefix for output files", required=False)
    parser.set_defaults(win=50, out="")
    return parser.parse_args()


def VAL_COUNTER(F1_seq_CHA, F2_seq_CHA):
    MUCH_num, Indel, Transver, Transition = 0, 0, 0, 0
    if F1_seq_CHA == F2_seq_CHA:
        MUCH_num += 1
    else:
        if F1_seq_CHA == "-" or F2_seq_CHA == "-":
            Indel +=1
        elif F1_seq_CHA == "A" and F2_seq_CHA == "C":
            Transver +=1
        elif F1_seq_CHA == "C" and F2_seq_CHA == "A":
            Transver +=1
        elif F1_seq_CHA == "G" and F2_seq_CHA == "T":
            Transver +=1
        elif F1_seq_CHA == "T" and F2_seq_CHA == "G":
            Transver +=1
        elif F1_seq_CHA == "A" and F2_seq_CHA == "G":
            Transition +=1
        elif F1_seq_CHA == "G" and F2_seq_CHA == "A":
            Transition +=1
        elif F1_seq_CHA == "C" and F2_seq_CHA == "T":
            Transition +=1
        elif F1_seq_CHA == "T" and F2_seq_CHA == "C":
            Transition +=1
    return MUCH_num, Indel, Transver, Transition
                    


nCr = {}
def cmb(n, r):
    if r == 0 or r == n: return 1
    if r == 1: return n
    if (n,r) in nCr: return nCr[(n,r)]
    nCr[(n,r)] = cmb(n-1,r) + cmb(n-1,r-1)
    return nCr[(n,r)]



if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    k_num = args.win
    plot_flag = args.plot
    out_prexi = args.out
    
    
    MUCH_num = 0
    UNMUCH_num = 0
    Indel =0 
    Transver =0
    Transition =0
    Already_done =[]

    label_total = ["Pair", "Consensus", "Indel", "Transversion", "Transition"]
    label = ["Pair", "start", "end", "Consensus", "Indel", "Transversion", "Transition"]
    labels = '\t'.join(label)
    
    
    df_total = pd.DataFrame(columns=label_total)
    df = pd.DataFrame(columns=label)
    
#full_length_analysis

    entries = [(fasta_in, fasta_in)]
    for index, (fasta_1, fasta_2) in enumerate(entries):
    ##fastaから配列を順番に読み込む
        with open(fasta_1, "rt") as fh1:
            for rec1 in SeqIO.parse(fh1, "fasta"):
                features=[] #featuresを初期化しておく(無い時があるので）
                length = len(rec1) #配列長を取得する
                F1_id = rec1.id #配列名を取得する
                F1_seq = rec1.seq.upper() #配列を大文字で取得する
                with open(fasta_2, "rt") as fh2:
                    for rec2 in SeqIO.parse(fh2, "fasta"):
                        features=[] #featuresを初期化しておく(無い時があるので）
                        F2_id = rec2.id #配列名を取得する
                        F2_seq = rec2.seq.upper() #配列を大文字で取得する
                        if F2_id != F1_id:
                            Pair = F1_id +  " vs " + F2_id
                            if Pair in Already_done:
                                for seq_num in range(length):
                                    tmp_MUCH_num, tmp_Indel, tmp_Transver, tmp_Transition  = VAL_COUNTER(F1_seq[seq_num], F2_seq[seq_num])
                                    MUCH_num += tmp_MUCH_num
                                    Indel    += tmp_Indel
                                    Transver += tmp_Transver
                                    Transition += tmp_Transition
                                sdf= pd.DataFrame(
                                                    {'Pair': Pair,
                                                    'Consensus': MUCH_num,
                                                    'Indel': Indel,
                                                    'Transversion': Transver,
                                                    'Transition': Transition}, index= [Pair] )
                                df_total = df_total.append(sdf)
                                MUCH_num = 0
                                Indel =0 
                                Transver =0
                                Transition =0
                            Already_done.append(F2_id+  " vs " + F1_id)
    print(df_total)
    df_total.to_csv(out_prexi + "whole_mutation_summary.txt", sep="\t")
    
    
    
#window_cut
    
    entries = [(fasta_in, fasta_in)]
    for index, (fasta_1, fasta_2) in enumerate(entries):
    ##fastaから配列を順番に読み込む
        with open(fasta_1, "rt") as fh1:
            for rec1 in SeqIO.parse(fh1, "fasta"):
                features=[] #featuresを初期化しておく(無い時があるので）
                length = len(rec1) #配列長を取得する
                F1_id = rec1.id #配列名を取得する
                F1_seq = rec1.seq.upper() #配列を大文字で取得する
                with open(fasta_2, "rt") as fh2:
                    for rec2 in SeqIO.parse(fh2, "fasta"):
                        features=[] #featuresを初期化しておく(無い時があるので）
                        F2_id = rec2.id #配列名を取得する
                        F2_seq = rec2.seq.upper() #配列を大文字で取得する
                        if F2_id != F1_id:
                            Pair = F1_id +  "vs" + F2_id
                            if Pair in Already_done:
                                for cut_length in range(0, length - k_num):
                                    F1_cut = F1_seq[cut_length:cut_length+k_num]
                                    F2_cut = F2_seq[cut_length:cut_length+k_num]
                                    for seq_num in range(k_num):
                                        tmp_MUCH_num, tmp_Indel, tmp_Transver, tmp_Transition  = VAL_COUNTER(F1_cut[seq_num], F2_cut[seq_num])
                                        MUCH_num += tmp_MUCH_num
                                        Indel    += tmp_Indel
                                        Transver += tmp_Transver
                                        Transition += tmp_Transition
                                        cut_length_mod = cut_length+k_num
                                    sdf_total= pd.DataFrame(
                                                        {'Pair': Pair,
                                                        'start': cut_length,
                                                        'end': cut_length_mod,
                                                        'Consensus': MUCH_num,
                                                        'Indel': Indel,
                                                        'Transversion': Transver,
                                                        'Transition': Transition}, index= [Pair + str(cut_length)] )
                                    df = df.append(sdf_total)
                                    MUCH_num = 0
                                    Indel =0 
                                    Transver =0
                                    Transition =0
                            Already_done.append(F2_id+  "vs" + F1_id)
                            
                            
                            

                                
    sns.set_style('whitegrid')
    
    df.start = df.start.astype(int)    
    df.end = df.end.astype(int)    
    df.Transversion = df.Transversion.astype(int)    
    df.Transition = df.Transition.astype(int)    
    df.Indel = df.Indel.astype(int)    
    df.Consensus = df.Consensus.astype(int)  
    
    df_mod =df[['start', 'Consensus','Indel','Transversion', 'Transition']]
    grouped = df_mod.groupby(['start'])
    groupde_mean = grouped.mean()
    
    groupde_mean.to_csv("average.txt", sep="\t")
    df.to_csv(out_prexi+"result_rawdata.txt", sep="\t")

    
    
    plt.figure()
    sns.countplot(x="Transition", data=df)
    plt.xlabel('Number of Transition par window')
    plt.savefig("hist_Transition.pdf")      
    if plot_flag == True:
        plt.show()
    
    plt.figure()
    sns.countplot(x="Transversion", data=df)
    plt.xlabel('Number of Transversion par window')
    plt.savefig(out_prexi+"hist_Transversion.pdf")      
    if plot_flag == True:
        plt.show()    
    
    plt.figure()
    sns.countplot(x="Indel", data=df)
    plt.xlabel('Number of Indel par window')
    plt.savefig(out_prexi+"hist_Indel.pdf")      
    if plot_flag == True:
        plt.show()

    plt.figure()
    g = sns.JointGrid(x="Transition", y ="Transversion", data=df)
    g.plot(sns.regplot, sns.distplot)
    plt.savefig("Transversion_vs_Transversion.pdf")
    plt.xlim(0,k_num)
    plt.ylim(0,k_num)
    if plot_flag == True:
        plt.show()

    plt.figure()
    sns.lineplot(x="start", y ="Transversion", hue = "Pair", data=df)
    plt.xlabel('start position of window')
    plt.savefig(out_prexi+"line_Transition.pdf")      
    if plot_flag == True:
        plt.show()    
    
    plt.figure()
    sns.lineplot(x="start", y ="Transition", hue = "Pair", data=df)
    plt.xlabel('start position of window')
    plt.savefig(out_prexi+"line_Transition.pdf")      
    if plot_flag == True:
        plt.show()    
    
    plt.figure()
    sns.lineplot(x="start", y ="Indel", hue = "Pair", data=df)
    plt.xlabel('start position of window')
    plt.savefig(out_prexi+"line_indel.pdf")      
    if plot_flag == True:
        plt.show()
    
    df_mod = df.melt(var_name = "mutation type", value_name = "mutation par window", id_vars = ["Pair", "start", "end"])

    
    plt.figure()
    sns.lineplot(x="start", y ="mutation par window", hue="mutation type", data=df_mod, ci="sd")
    plt.xlabel('start position of window')
    plt.ylabel('mutation par window (Mesh = SD)')
    plt.savefig(out_prexi+"line_all.pdf")      
    if plot_flag == True:
        plt.show()