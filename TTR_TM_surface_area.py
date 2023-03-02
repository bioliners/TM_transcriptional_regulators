# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 08:28:33 2021

@author: LukeD
"""

import os
import csv

os.chdir(os.path.dirname(os.path.abspath(__file__)))

end_of_prot=0
i=0
prot_seq=[]
prot_ID=[]
total_protein_list=[]
total_num_proteins=0
line_counter=0
TM_len=[]
TM_sur_area=0

'''amino acid surface area calculations here were preformed as previously described (1)
1. Yuan, Z., Zhang, F., Davis, M.J., Bodén, M. and Teasdale, R.D., 2006. Predicting the solvent accessibility of transmembrane residues from protein sequence. Journal of proteome research, 5(5), pp.1063-1070.'''
sur_area_ind ={"A":116.4, "R":249.26, "N":168.87, "D":155.37, "C":141.48, "Q":189.17, 
                 "E":187.16, "G":83.91, "H":198.51, "I":189.95,"L":197.99, "K":207.49, 
                 "M":210.55, "F":223.29, "P":144.8, "S":125.68, "T":148.06, "W":265.42, 
                 "Y":238.3, "V":162.24, "\n":0, 'X':0}

with open("TTR_domains_output.csv", "w", newline='') as TTR_TM_analysis:
    wr = csv.writer(TTR_TM_analysis, quoting=csv.QUOTE_ALL)
    genome_file= open('TTR_domains.fa','r')
    def genome_file_len(genome_file):
        with open('TTR_domains.fa','r') as f:
            for i, l in enumerate(f):
                pass
            return i+1
    #print('number of lines: ', genome_file_len(genome_file))
    
    for x in range(0,genome_file_len(genome_file)):
        genome_line=genome_file.readline()
        genome_line.strip('\n')
        if '>' in genome_line:
            prot_ID.append(genome_line[1:])
            line_counter+=1
            
        if '>' not in genome_line:
            line_counter=0
            TM_len=[]
            line_len=0
            line_len= len(genome_line)-1
            TM_len.append(line_len)
            TM_sur_area=0
            TM_sur_area_list=[]
            TM_sur_list=[]
            A_list=[]
            C_list=[]
            D_list=[]
            E_list=[]
            F_list=[]
            G_list=[]
            H_list=[]
            I_list=[]
            K_list=[]
            L_list=[]
            M_list=[]
            N_list=[]
            P_list=[]
            Q_list=[]
            R_list=[]
            S_list=[]
            T_list=[]
            V_list=[]
            W_list=[]
            Y_list=[]
            A_counter=0
            C_counter=0
            D_counter=0
            E_counter=0
            F_counter=0
            G_counter=0
            H_counter=0
            I_counter=0
            K_counter=0
            L_counter=0
            M_counter=0
            N_counter=0
            P_counter=0
            Q_counter=0
            R_counter=0
            S_counter=0
            T_counter=0
            V_counter=0
            W_counter=0
            Y_counter=0
            for x in genome_line:
                TM_sur_area_list.append(sur_area_ind[x])
                if x == 'A':
                    A_counter+=1
                if x == 'C':
                    C_counter+=1  
                if x == 'D':
                    D_counter+=1
                if x == 'E':
                    E_counter+=1
                if x == 'F':
                    F_counter+=1
                if x == 'G':
                    G_counter+=1
                if x == 'H':
                    H_counter+=1   
                if x == 'I':
                    I_counter+=1    
                if x == 'K':
                    K_counter+=1  
                if x == 'L':
                    L_counter+=1  
                if x == 'M':
                    M_counter+=1
                if x == 'N':
                    N_counter+=1
                if x == 'P':
                    P_counter+=1
                if x == 'Q':
                    Q_counter+=1
                if x == 'R':
                    R_counter+=1
                if x == 'S':
                    S_counter+=1
                if x == 'T':
                    T_counter+=1
                if x == 'V':
                    V_counter+=1
                if x == 'W':
                    W_counter+=1
                if x == 'Y':
                    Y_counter+=1
            A_list.append(A_counter)
            C_list.append(C_counter)
            D_list.append(D_counter)
            E_list.append(E_counter)
            F_list.append(F_counter)
            G_list.append(G_counter)
            H_list.append(H_counter)
            I_list.append(I_counter)
            K_list.append(K_counter)
            L_list.append(L_counter)
            M_list.append(M_counter)
            N_list.append(N_counter)
            P_list.append(P_counter)
            Q_list.append(Q_counter)
            R_list.append(R_counter)
            S_list.append(S_counter)
            T_list.append(T_counter)
            V_list.append(V_counter)
            W_list.append(W_counter)
            Y_list.append(Y_counter)
            TM_sur_area= sum(TM_sur_area_list)
            TM_sur_list.append(TM_sur_area)
            prot_ID=prot_ID[0:] + TM_len + TM_sur_list + A_list + C_list + D_list + E_list + F_list + G_list + H_list + I_list + K_list + L_list + M_list + N_list + P_list + Q_list + R_list + S_list + T_list + V_list + W_list + Y_list
            wr.writerow(prot_ID)
            print(prot_ID)
            #print(prot_ID)
            prot_ID=[]
genome_file.close()
       
