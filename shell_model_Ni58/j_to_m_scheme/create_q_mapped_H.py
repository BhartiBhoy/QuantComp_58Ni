from .j_to_m_Hamiltonian import Hamiltonian_M_scheme

from fractions import Fraction
import numpy as np
from math import sqrt
import math
from itertools import combinations
import os 

def Read_model_space(model_filename):
    current_dir=os.getcwd()
    path = 'shell_model_Ni58/model_space'
    full_path=os.path.join(os.path.join(current_dir,path),model_filename)

    with open(full_path, 'r') as f:
        lines = f.readlines()
        l1 = lines[0].split('=')[1].replace('[','').replace(']','').split(',')
        # print(l1)
        levels = {}
        for l in l1:
            data = l.split(':')
            # print(data)
            levels[int(data[0])] = (data[1].replace('\'', '')).replace('\n', '')
        print('levels = {}'.format(levels))


        l2 = lines[1].split('=')[1].replace('[', '').replace(']', '').split(',')
        # print(l2)
        energies = {}
        for l in l2:
            data = l.split(':')
            energies[int(data[0])] = float(data[1])
        print('energes = {}'.format(energies))

    return levels,energies


def orbital_to_l(orb):
        orb_to_l_dic = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'j': 7}
        return orb_to_l_dic[orb]

def All_Ms(levels):
        Ms = {}
        for key, val in levels.items():
            val = val.split(';')[0]
            j = float(Fraction(val[2:]))
            Ms[key] = []
            for m in np.arange(j, -j-1, -1):
                Ms[key].append(m)
            # print(Ms[key])

        k1 = -1
        Ms_all=[]
        for key1, m1s in Ms.items():
            k1 += 1
            # print(k1)
            for i, m1 in enumerate(m1s):
                for key2, m2s in list(Ms.items())[k1:]:
                    if key1 == key2:
                        j = i+1
                        for m2 in m2s[i+1:]:
                            Ms_all.append(m1+m2)
                
                            j += 1
                    else:
                        for j, m2 in enumerate(m2s):
                            Ms_all.append(m1+m2)
                      
        M_vals = -np.sort(-np.array(list(set(Ms_all))))
        return M_vals



def Read_interaction_J(interaction_filename,levels):
    Correct_Qbit_map = {}
    for key, val in levels.items():

        Correct_Qbit_map[val] = str(key)
    
    current_dir=os.getcwd()
    path = 'shell_model_Ni58/interactions'
    full_path=os.path.join(os.path.join(current_dir,path),interaction_filename)
    with open(full_path,'r') as f:
        
        lines = f.readlines()
        line1=lines[0]
        data = line1.split(',')
   
        interaction_levels={}
        for d in data:
            key = d.split('=')[0].strip()
            val = d.split('=')[1].strip()
            interaction_levels[key] = val
        lines = lines[9:]
        data = lines[0].split()
        D_j_scheme={}
        for lineno,line in enumerate(lines[1:]):
            data = line.split()
            # print(data)
            if len(data)>1:
                # print(lineno)
                a,b,c,d,J,T,val = data[0],data[1],data[2],data[3],data[4],data[5],\
                                float(data[6])
                a,b,c,d = Correct_Qbit_map[interaction_levels[a]],Correct_Qbit_map[interaction_levels[b]],Correct_Qbit_map[interaction_levels[c]],Correct_Qbit_map[interaction_levels[d]]


            key = a+' '+b+' '+c+' '+d+' '+J+' '+T
            D_j_scheme[key] = val
        
        for lineno,line in enumerate(lines[1:]):
            data = line.split()
            if len(data)>1:
                a,b,c,d,J,T,val = data[0],data[1],data[2],data[3],data[4],data[5],\
                                float(data[6])
                a,b,c,d = Correct_Qbit_map[interaction_levels[a]],Correct_Qbit_map[interaction_levels[b]],Correct_Qbit_map[interaction_levels[c]],Correct_Qbit_map[interaction_levels[d]]

            keys = [a+' '+b+' '+d+' '+c+' '+J+' '+T,b+' '+a+' '+c+' '+d+' '+J+' '+T,b+' '+a+' '+d+' '+c+' '+J+' '+T,
                    c+' '+d+' '+a+' '+b+' '+J+' '+T,d+' '+c+' '+a+' '+b+' '+J+' '+T,c+' '+d+' '+b+' '+a+' '+J+' '+T,d+' '+c+' '+b+' '+a+' '+J+' '+T]

            for k in keys:
                if k not in D_j_scheme.keys():
                    # print(k,val)
                    D_j_scheme[k] = val


        return D_j_scheme

def Create_Qmap(levels):
    Qubit_map = {}
    qbit_no = 0
    for orbital in levels.values():
        n = int(orbital[0])
        l = orbital_to_l(orbital[1])
        j = float(Fraction(orbital.split(';')[0][2:]))
        nuc = orbital.split(';')[1]
        Jz = j
        while Jz>0:

            key = str(n)+','+str(l)+','+str(j)+','+str(-Jz)+','+nuc
            Qubit_map[key] = qbit_no
            key = str(n)+','+str(l)+','+str(j)+','+str(Jz)+','+nuc
            Qubit_map[key] = qbit_no+1
            qbit_no = qbit_no+2
            Jz = Jz-1

    return Qubit_map



def Write_Qmapped_M_Hamiltonian(levels,energies,M_vals,Interaction_Dict,Qubit_map,output_filename):

    current_dir=os.getcwd()
    path = 'shell_model_Ni58/q_mapped_H_Output'
    full_path=os.path.join(os.path.join(current_dir,path),output_filename)
    Nuc_Tz = {0.5:'n',-0.5:'p'}
    ''' printing output in the output file'''
    with open(full_path, 'w') as f:
        f.write('Qubit Map of the M-scheme states \n')
        f.write('No. of Qubits = {} \n'.format(len(Qubit_map.keys())))
        # f.write('________________________________________________\n')
        f.write('Qubit\tn\tl\tj\t  jz \n')
        print('Qubit    n       l        j      jz \n')
        
        for key,val in Qubit_map.items():
            nlj_jz = key.split(',')
            print('  ',val,'\t',nlj_jz[0],'\t',nlj_jz[1],'\t',Fraction(float(nlj_jz[2])),'\t',Fraction(float(nlj_jz[3])))
            f.write('\t%s\t%s\t%s\t%s\t  %s \n'%(val,nlj_jz[0],nlj_jz[1],Fraction(float(nlj_jz[2])),Fraction(float(nlj_jz[3]))))
        # All_M_Bases,All_V,All_KE=[],[],[]
        f.write('\n')
        f.write('________________________________________________\n')
        f.write('Qubit mapped matrix elements of V All M \n')
        print('________________________________________________\n')
        # print('Qubit mapped matrix elements of V All M \n')

        countM=0
        for M in M_vals:
            # print(M)

            '''Getting the m-scheme basis and the Hamiltonian'''
            Base_WFs,V,KE = Hamiltonian_M_scheme(levels,energies,M,Interaction_Dict,Qubit_map)
            for i,base in enumerate(Base_WFs):
                kE_sp=energies[base[5]]+energies[base[11]]
                V[i,i] += kE_sp
            Qbit_bases=[]
            for base in Base_WFs:
                n1,l1,j1,jz1,tz1,n2,l2,j2,jz2,tz2 = base[0],base[1],base[2],base[3],base[4],base[6],base[7],base[8],base[9],base[10]
                nuc1 = Nuc_Tz[tz1]
                nuc2 = Nuc_Tz[tz2]
                key1 = str(n1)+','+str(l1)+','+str(j1)+','+str(jz1)+','+str(nuc1)
                key2 = str(n2)+','+str(l2)+','+str(j2)+','+str(jz2)+','+str(nuc2)
                Qbit1 = Qubit_map[key1]
                Qbit2 = Qubit_map[key2]
                Qbit_bases.append([Qbit1,Qbit2])
            

            for i,qbit_base1 in enumerate(Qbit_bases):
                for j,qbit_base2 in enumerate(Qbit_bases):
                    if abs(V[i,j])>0:
                        # print(qbit_base1[0],qbit_base1[1],qbit_base2[0],qbit_base2[1],'  ',V[i,j])
                        f.write('   %d\t%d\t%d\t%d\t\t%f \n'%(qbit_base1[0],qbit_base1[1],qbit_base2[0],qbit_base2[1],V[i,j]))
                        # f.write('   %d\t%d\t%d\t%d\t\t%f \n'%(qbit_map_trans[qbit_base1[0]],qbit_map_trans[qbit_base1[1]],qbit_map_trans[qbit_base2[0]],qbit_map_trans[qbit_base2[1]],V[i,j]))
                        countM+=1
        f.write('End V and Total Matrix Elements = {}'.format(countM))

        
        
        f.write('\n')
        f.write('________________________________________________\n')
        f.write('Qubit mapped matrix elements of KE all M \n')
        # print('________________________________________________\n')
        # print('Qubit mapped matrix elements of KE for All M \n')
        for M in M_vals:
            # print(M)

            '''Getting the m-scheme basis and the Hamiltonian'''
            Base_WFs,V,KE = Hamiltonian_M_scheme(levels,energies,M,Interaction_Dict,Qubit_map)

            Qbit_bases=[]
            for base in Base_WFs:
                n1,l1,j1,jz1,tz1,n2,l2,j2,jz2,tz2 = base[0],base[1],base[2],base[3],base[4],base[6],base[7],base[8],base[9],base[10]
          
                nuc1 = Nuc_Tz[tz1]
                nuc2 = Nuc_Tz[tz2]
                key1 = str(n1)+','+str(l1)+','+str(j1)+','+str(jz1)+','+str(nuc1)
                key2 = str(n2)+','+str(l2)+','+str(j2)+','+str(jz2)+','+str(nuc2)
                Qbit1 = Qubit_map[key1]
                Qbit2 = Qubit_map[key2]
                Qbit_bases.append([Qbit1,Qbit2])
            for i,qbit_base1 in enumerate(Qbit_bases):
                # print(qbit_base1[0],qbit_base1[1],qbit_base1[0],qbit_base1[1],'  ',KE[i,i])
                f.write('   %d\t%d\t%d\t%d\t\t%f \n'%(qbit_base1[0],qbit_base1[1],qbit_base1[0],qbit_base1[1],KE[i,i]))
        f.write('End KE')