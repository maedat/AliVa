#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : AliVa
# File Name       : aliva.py
# Version       : 0.0.1
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
from BCBio import GFF





def GET_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="File path to a genome sequence file", required=True)
    parser.add_argument('-k','--kmer', help="k-mear number", required=False)
    parser.set_defaults(kmer=50)
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
                    






if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    k_num = args.kmer
    
    
    MUCH_num = 0
    UNMUCH_num = 0
    Indel =0 
    Transver =0
    Transition =0

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
                            for cut_length in range(0, length - k_num):
                                F1_cut = F1_seq[cut_length:cut_length+k_num]
                                F2_cut = F2_seq[cut_length:cut_length+k_num]
                                for seq_num in range(k_num):
                                    tmp_MUCH_num, tmp_Indel, tmp_Transver, tmp_Transition  = VAL_COUNTER(F1_cut[seq_num], F2_cut[seq_num])
                                    MUCH_num += tmp_MUCH_num
                                    Indel    += tmp_Indel
                                    Transver += tmp_Transver
                                    Transition += tmp_Transition
                                OUT = [Pair, str(cut_length), str(cut_length+k_num), str(MUCH_num), str(Indel), str(Transver), str(Transition)]
                                OUTs = '\t'.join(OUT)
                                print(OUTs)
                                MUCH_num = 0
                                Indel =0 
                                Transver =0
                                Transition =0
