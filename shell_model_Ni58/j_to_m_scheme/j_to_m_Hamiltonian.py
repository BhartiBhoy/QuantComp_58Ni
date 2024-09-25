
from fractions import Fraction
import numpy as np
import math

def Hamiltonian_M_scheme(levels,energies,M,Interaction_Dict,Qubit_map):
    Nuc_Tz = {0.5:'n',-0.5:'p'}
    Tz_dict = {'p':-0.5,'n':0.5}
    def Gen_basis(levels, M):
        ''' Generate m-scheme basis and projections of js in an array Ms from the dicionery levels and the total M'''
        Ms = {}
        for key, val in levels.items():
    
            j = float(Fraction(val.split(';')[0][2:]))
            Ms[key] = []
            for m in np.arange(j, -j-1, -1):
                Ms[key].append(m)
     

        k1 = -1
        basis = []
        for key1, m1s in Ms.items():
            k1 += 1
            for i, m1 in enumerate(m1s):
                for key2, m2s in list(Ms.items())[k1:]:
                    if key1 == key2:
                        j = i+1
                        for m2 in m2s[i+1:]:
                            if m1+m2 == M:
                                basis.append([key1, i, key2, j])
                            j += 1
                    else:
                        for j, m2 in enumerate(m2s):
                            if m1+m2 == M:
                                basis.append([key1, i, key2, j])

        return Ms, basis  

    def orbital_to_l(orb):
        ''' Obtain the  l value corresponding to the symbols (s,p,d,f,..) '''
        orb_to_l_dic = {'s': 0, 'p': 1, 'd': 2,
                        'f': 3, 'g': 4, 'h': 5, 'i': 6, 'j': 7}
     
        return orb_to_l_dic[orb]


    def notation_input(string_a):
        '''Coverts notation like 1s1/2 to respective n,l,j values'''
        from fractions import Fraction
        # string_a = '0s1/2'
        na = int(string_a[0])
        la = orbital_to_l(string_a[1])
        ja = float(Fraction(string_a.split(';')[0][2:]))
        tz = Tz_dict[string_a.split(';')[1]]
   
        return na, la, ja,tz
    
    def NLJM_values(i, Ms,base):
        '''Obtain n,l,j,m values using the array basis'''
        level = levels[base[i]]
        n, l, j,tz = notation_input(level)
        m = Ms[base[i]][base[i+1]]
      
        return n, l, j, m,tz
    
    Ms,basis = Gen_basis(levels,M)
 
    # Generating the full m-scheme basis [[na,la,ja,ma,a, nb,lb,jb,mb,b],[..],...]
    M_scheme_basis = []
    for i, base in enumerate(basis):
        na, la, ja, ma,tza = NLJM_values(0,Ms, base)
        nb, lb, jb, mb,tzb = NLJM_values(2,Ms, base)
        temp = [na, la, ja, ma,tza, base[0], nb, lb, jb, mb,tzb, base[2]]
        if ma+mb == M:
            M_scheme_basis.append(temp)
        else:
            print('Invalid Input in basis number {}'.format(i+1))
    
    
    def fact(n):
      
        n = int(n)
        return (math.factorial(n))

    def delta(a, b):
        if a == b:
            return 1
        else:
            return 0
    
    def cgc(j1, j2, m1, m2, J):
        M = m1+m2
        if (abs(m1+m2) > J):
            return 0

        term1a = (2*J+1)*fact(J+j1-j2)*fact(J-j1+j2)*fact(j1+j2-J)
        term1b = fact(j1+j2+J+1)
        term1 = math.sqrt(term1a/term1b)

        term2a = fact(J+M)*fact(J-M)*fact(j1-m1) * \
            fact(j1+m1)*fact(j2-m2)*fact(j2+m2)
        term2 = math.sqrt(term2a)

        ''' Some conditions on integer k for which every factorial should be nonnegative '''

        A = j1+j2-J
        B = j1-m1
        C = j2+m2
        D = J-j2+m1
        E = J-j1-m2

        minik = np.minimum(D, E)
        maxik = np.min([A, B, C])

        if minik < 0.0:
            minik = abs(minik)
        else:
            minik = 0.0

        term3 = 0.0
        for k in range(int(minik), int(maxik+1)):
            term3a = (-1)**k
            term3b = fact(k)*fact(j1+j2-J-k)*fact(j1-m1-k) * \
                    fact(j2+m2-k)*fact(J-j2+m1+k)*fact(J-j1-m2+k)
            term3c = term3a/term3b
            term3 = term3 + term3c

        answer = term1*term2*term3
        return answer
    
    def interaction(a,b,c,d,J,T):
        key = str(int(a))+' '+str(int(b))+' '+str(int(c))+' '+str(int(d))+' '+str(int(J))+' '+str(int(T))
        if key in Interaction_Dict.keys():
            return Interaction_Dict[key]
        else:
            return 0



    def T_e(a, b, c, d, J, T, e_a, e_b):
        '''Calculates the KE matrix element in JT scheme'''
        Const1 = delta(a, c)*delta(b, d)
        Const2 = e_a+e_b
        Const = Const1*Const2
        answer = Const*(1-delta(a, b)*((-1)**(J+T)))/(1+delta(a, b))
        return answer

    def CG(j1, m1, j2, m2, J, M):
        ans = cgc(j1, j2, m1, m2, J)
        return ans
    
    
    def M_scheme(M_s_b):
        '''Calculate the Hamiltonian Matrix in m-scheme, the m-scheme basis array is the input'''

        # Hamiltonian of dimension N*N where N is the basis size
        H_mscheme = np.zeros((len(M_s_b), len(M_s_b)), dtype=float)
        KE_mscheme = np.zeros((len(M_s_b), len(M_s_b)), dtype=float)
        allowed={}
        for row, m_bas_1 in enumerate(M_s_b):
            na, la, ja, ma, tza,a, nb, lb, jb, mb,tzb, b = m_bas_1
            nuca = Nuc_Tz[tza]
            nucb = Nuc_Tz[tzb]
            keya = str(na)+','+str(la)+','+str(ja)+','+str(ma)+','+str(nuca)
            keyb = str(nb)+','+str(lb)+','+str(jb)+','+str(mb)+','+str(nucb)
            Qbita = Qubit_map[keya]
            Qbitb = Qubit_map[keyb]

            if ja == jb and a == b:             # checks whether two levels are same and two particles have same J number
                Js1 = [i for i in np.arange(0, ja+jb+1, 2)]
                Ts1=[1]
            else:
                Js1 = [i for i in np.arange(abs(ja-jb), abs(ja+jb)+1, 1)]
                if tza==tzb:
                    Ts1=[1]
                else:
                    Ts1=[i for i in np.arange(abs(tza-tzb), abs(tza+tzb)+1, 1)]
            
            
            for col, m_bas_2 in enumerate(M_s_b):
                nc, lc, jc, mc, tzc, c, nd, ld, jd, md, tzd, d = m_bas_2
                nuc_c = Nuc_Tz[tzc]
                nucd = Nuc_Tz[tzd]
                keyc = str(nc)+','+str(lc)+','+str(jc)+','+str(mc)+','+str(nuc_c)
                keyd = str(nd)+','+str(ld)+','+str(jd)+','+str(md)+','+str(nucd)
                Qbitc = Qubit_map[keyc]
                Qbitd = Qubit_map[keyd]
                # qno = [na, la, ja, nb, lb, jb, nc, lc, jc, nd, ld, jd]
                if jc == jd and c == d:
                    Js2 = [i for i in np.arange(0, jc+jd+1, 2)]
                    Ts2=[1]
                else:
                    Js2 = [i for i in np.arange(abs(jc-jd), abs(jc+jd)+1, 1)]
                    if tzc==tzd:
                        Ts2=[1]
                    else:
                        Ts2=[i for i in np.arange(abs(tzc-tzd), abs(tzc+tzd)+1, 1)]
           
                H_mscheme[row, col] = 0.0
                KE_mscheme[row, col] = 0.0
               
                
                for J1 in Js1:
                    for J2 in Js2:
                        if J1 == J2:
                            for T1 in Ts1:
                                for T2 in Ts2:
                                    if T1==T2 and (tza+tzb)==(tzc+tzd):
      
                                        # gets sp energies depening on a and b
                                        e_a, e_b = energies[a], energies[b]
                                        # Kinetic Energy Matrix element with T=1 to reduce into the formula for J only
                                        KE = T_e(a, b, c, d, J1, 1, e_a, e_b)
                                        N_ab = np.sqrt((1-delta(a, b)*(-1)**(J1+T1)))/(1+delta(a, b))
                                        N_cd = np.sqrt((1-delta(c, d)*(-1)**(J2+T2)))/(1+delta(c, d))
                                        N = N_ab*N_cd
                                        C_G_coef = CG(ja, ma, jb, mb, J1, M)*CG(jc, mc, jd, md, J2, M)*CG(0.5,tza,0.5,tzb,T1,tza+tzb)*CG(0.5,tzc,0.5,tzd,T2,tzc+tzd)
                          
                                        if interaction(a,b,c,d,J1,T1) !=0:
                                            H_mscheme[row, col] += C_G_coef *(interaction(a,b,c,d,J1,T1))/N
                                        else:
                                            continue
                                        KE_mscheme[row, col] += C_G_coef*(KE)/N

        return H_mscheme,KE_mscheme


    V_M,KE_M = M_scheme(M_scheme_basis)   # calculating the Hamiltonian
    
    return M_scheme_basis,V_M,KE_M
