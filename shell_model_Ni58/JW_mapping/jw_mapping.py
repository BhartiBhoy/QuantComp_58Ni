import numpy as np
import qiskit_nature
from qiskit.primitives import Estimator, BackendEstimator
from qiskit.algorithms.minimum_eigensolvers import VQE
from qiskit.algorithms.optimizers import COBYLA, SLSQP,SPSA,GradientDescent
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.utils import algorithm_globals
from qiskit import QuantumCircuit
from qiskit.circuit import Gate
from qiskit import Aer


import matplotlib.pyplot as plt

from fractions import Fraction
import os 


qiskit_nature.settings.use_pauli_sum_op = False


def Read_m_scheme_H(m_scheme_H_filename,levels,energies):
    current_dir=os.getcwd()
    path = 'shell_model_Ni58/q_mapped_H_output'
    full_path=os.path.join(os.path.join(current_dir,path),m_scheme_H_filename)
    # print(full_path)

        
    with open(full_path,'r') as f:
        lines = f.readlines()
        N_Qbits = int(lines[1].split('=')[1])
        # print(lines)
        two_body_string = 'Qubit mapped matrix elements of V All M'
        two_body_string_end = 'End V'

        One_body_string = 'Qubit mapped matrix elements of KE all M'
        One_body_string_end = 'End KE'

        for i,line in enumerate(lines):
            if two_body_string in line:
                # print(i,line)
                start_V = i
            if two_body_string_end in line:
                # print(i,line)
                end_V = i
            if One_body_string in line:
                # print(i,line)
                start_KE = i
            if One_body_string_end in line:
                # print(i,line)
                end_KE = i
        
        Hamiltonian_Coeff = {}
        for line in lines[start_V+1:end_V]:
            # print(line)
            data = line.split()
            key = '+_{} +_{} -_{} -_{}'.format(data[0],data[1],data[3],data[2]) #0=i,1=j,2=k,3=l
            value = float(data[-1])
            # print(key,value)
            Hamiltonian_Coeff[key] = value#*0.5        #multiply by 0.5
        
        for line in lines[start_KE+1:end_KE]:
            data = line.split()
            key = '+_{} -_{}'.format(data[0],data[1])
            value = float(data[-1])
            # print(key,value)
            # Hamiltonian_Coeff[key] = value


    ene=[]
    for key, val in levels.items():
        j = float(Fraction(val.split(';')[0][2:]))
        ene = ene+[energies[key]]*int(2*j+1)
    # print(ene)

    for i in range(N_Qbits):
        key = '+_{} -_{}'.format(i,i)
        value = ene[i]
        # print(key,value)
        # Hamiltonian_Coeff[key] = value


    # print(end_V-start_V-1)


    return Hamiltonian_Coeff,N_Qbits


def Create_JW_mapping(H_coeff,N_Qbits):
    
    Hamiltonian_op = FermionicOp(
        H_coeff,
        num_spin_orbitals=N_Qbits,
    )
    
    mapper = JordanWignerMapper()
    qubit_hamiltonian = mapper.map(Hamiltonian_op)
    if not isinstance(qubit_hamiltonian, SparsePauliOp):
            qubit_hamiltonian = qubit_hamiltonian.primitive
    

    return qubit_hamiltonian


            


def Single_excitation(phi):
    circ = QuantumCircuit(2,name='Single Excitation')
    circ.cx(0,1)
    circ.ry(phi/2,0)
    circ.cx(1,0)
    circ.ry(-phi/2,0)
    circ.cx(1,0)
    circ.cx(0,1)
    return circ


def Double_excitation(theta):
    circ = QuantumCircuit(4,name='Double Excitation')
    circ.cx(2,3)
    circ.cx(0,2)
    circ.h([0,3])
    circ.cx(0,1)
    circ.cx(2,3)
    circ.ry(-theta/8,0)
    circ.ry(theta/8,1)
    circ.cx(0,3)
    circ.h(3)
    circ.cx(3,1)
    circ.ry(-theta/8,0)
    circ.ry(theta/8,1)
    circ.cx(2,1)
    circ.cx(2,0)
    circ.ry(theta/8,0)
    circ.ry(-theta/8,1)
    circ.cx(3,1)
    circ.h(3)
    circ.cx(0,3)
    circ.ry(theta/8,0)
    circ.ry(-theta/8,1)
    circ.cx(0,1)
    circ.cx(2,0)
    circ.h([0,3])
    circ.cx(0,2)
    circ.cx(2,3)
    return circ


def Ensatz_ground(Qbits,InitialState,Thetas):
    circ = QuantumCircuit(Qbits,name='Ensatz')
    circ.x(InitialState)
    Nmax_dble = Qbits//2 - 1
    Initial_qbit=0
    for i in range(Nmax_dble):
        doubleExci = Double_excitation(Thetas[i])
        circ.append(doubleExci,list(range(Initial_qbit,Initial_qbit+4)))
        Initial_qbit += 2
    return circ


def Ensatz_first(Qbits,InitialState,Phis,Ground,Excite,double=True):
    # print(Ground,Excite)
    # print(len(Ground),len(Excite))
    circ = QuantumCircuit(Qbits,name='Exci_Ensatz')
    circ.x(InitialState)
    # for i in range(2)
    singleExcite = Single_excitation(Phis[0]) # 12 qubit
    # circ.append(singleExcite,[3,5])
    circ.append(singleExcite,[3,9])

    singleExcite = Single_excitation(Phis[1]) # 12 qubit
    circ.append(singleExcite,[3,11])

    singleExcite = Single_excitation(Phis[2]) # 12 qubit
    circ.append(singleExcite,[1,7])


    if double:
        for i in range(len(Ground)):
            doubleExci = Double_excitation(Phis[i+2])
            circ.append(doubleExci,Ground[i]+Excite[i])
    
    return circ


def Ensatz_second(Qbits,InitialState,Phis,Ground,Excite,double=True):
    print(Ground,Excite)
    print(len(Ground),len(Excite))
    circ = QuantumCircuit(Qbits,name='Exci_Ensatz')
    circ.x(InitialState)
    # for i in range(2)
    singleExcite = Single_excitation(Phis[0]) # 12 qubit
    circ.append(singleExcite,[1,7])
    if double:
        for i in range(len(Ground)):
            doubleExci = Double_excitation(Phis[i+2])
            circ.append(doubleExci,Ground[i]+Excite[i])
    
    return circ


def Ansatz(N_Qbits,state,output_fig_path):

    current_dir=os.getcwd()
    
    full_path=os.path.join(current_dir,output_fig_path)

    if state==0:
        Parameters=[]
        for i in range(N_Qbits//2-1):
            theta = Parameter('theta_{}'.format(i))         # for different theta at every excitation
            Parameters.append(theta)
    
        gs_ini = [0,1]
        ansatz = Ensatz_ground(N_Qbits,gs_ini,Parameters)
        ansatz.draw('mpl').savefig(full_path+'/Circuit_for_Ensatz_ground.png')
        print(f'Circuit Figure Saved in {full_path} directory')
    elif state==1:        
        Parameters=[]
        for i in range(23):
            theta = Parameter('theta_{}'.format(i))         # for different theta at every excitation
            Parameters.append(theta)
        gs_ini = [1,3]
        no_of_ops = 5
        db_exc = [[5,8],[5,10],[2,5],[7,9],[7,11]]
        ansatz = Ensatz_first(N_Qbits,gs_ini,Parameters,Ground=[gs_ini]*no_of_ops,Excite=db_exc)  # 12 qubit
        ansatz.draw('mpl').savefig(full_path+'/Circuit_for_Ensatz_first_excited.png')
        print(f'Circuit Figure Saved in {full_path} directory')
    elif state==2:        
        Parameters=[]
        for i in range(23):
            theta = Parameter('theta_{}'.format(i))         # for different theta at every excitation
            Parameters.append(theta)
        gs_ini = [1,5]
        no_of_ops = 0
        db_exc = []
        ansatz = Ensatz_second(N_Qbits,gs_ini,Parameters,Ground=[gs_ini]*no_of_ops,Excite=db_exc)  # 12 qubit
        ansatz.draw('mpl').savefig(full_path+'/Circuit_for_Ensatz_second_excited.png')
        print(f'Circuit Figure Saved in {full_path} directory')

    return ansatz



def VQE_run(ansatz,qubit_hamiltonian):
    
    seed = 500
    algorithm_globals.random_seed = seed
    # for i in range(runs)
    counts = []
    values = []
    params = []
    deviation = []
    

    def callback(eval_count, parameters, mean, std):
        # Overwrites the same line when printing
        print("Evaluation: {}, Energy: {}, Std: {}".format(eval_count, mean, std))
        # clear_output(wait=True)
        counts.append(eval_count)
        values.append(mean)
        params.append(parameters)
        deviation.append(std)
       


    Iterations=1000
    # backend = Aer.get_backend('statevector_simulator')
    # backend = Aer.get_backend('qasm_simulator')

    Optimizer_name = 'COBYLA'
    # Optimizer_name = 'SPSA'
    # Optimizer_name = 'SLSQP'
    # Optimizer_name = 'GradientDescent'

    if Optimizer_name == 'COBYLA':
        optimizer = COBYLA(maxiter=Iterations)
    elif Optimizer_name == 'SPSA':
        optimizer = SPSA(maxiter=Iterations)
    elif Optimizer_name == 'SLSQP':
        optimizer = SLSQP(maxiter=Iterations)
    elif Optimizer_name == 'GradientDescent':
        optimizer = GradientDescent(maxiter=Iterations, learning_rate=0.1)

    vqe = VQE(
        Estimator(),
        # BackendEstimator(backend=backend,options={'shots':8000}),
        ansatz=ansatz,
        optimizer=optimizer,
        callback=callback,
    )

        
    vqe_result = vqe.compute_minimum_eigenvalue(qubit_hamiltonian)
    errors =[]
    for d in deviation:
        if 'variance' in d.keys():
            errors.append(np.sqrt(d['variance'])/2)
    # if len(errors)>1:
    #     for c,v,e in zip(counts,values,errors):
    #         print('{} \t {} \t {} \n'.format(c,v,e))
    # else:
    #     for c,v in zip(counts,values):
    #         print('{} \t {} \n'.format(c,v))