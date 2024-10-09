from shell_model_Ni58.j_to_m_scheme.create_q_mapped_H import Read_model_space,All_Ms,Read_interaction_J,Create_Qmap
from shell_model_Ni58.j_to_m_scheme.create_q_mapped_H import Write_Qmapped_M_Hamiltonian
from shell_model_Ni58.JW_mapping.jw_mapping import Read_m_scheme_H,Create_JW_mapping,Ansatz,VQE_run

import numpy as np



model_filename='model_space.in'
interaction_filename = 'jun45_pf.int'
output_file_for_qmapped_H='Qbit_mapped_all_M_JUN45.out'

state=0 # ground state
state=1 # 1st excited state
state=2 # 2nd excited state

levels,energies=Read_model_space(model_filename)

M_vals=All_Ms(levels)

print(M_vals)

Interaction_Dict=Read_interaction_J(interaction_filename,levels)

# print(Interaction_Dict)
Qubit_map = Create_Qmap(levels)

# print(Qubit_map)
Write_Qmapped_M_Hamiltonian(levels,energies,M_vals,Interaction_Dict,Qubit_map,output_file_for_qmapped_H)

Hamiltonian_Coeff,N_qbits = Read_m_scheme_H(m_scheme_H_filename=output_file_for_qmapped_H,levels=levels,energies=energies)

qubit_hamiltonian = Create_JW_mapping(Hamiltonian_Coeff,N_qbits)

ansatz = Ansatz(N_Qbits=N_qbits,state=state,output_fig_path='Figures')

VQE_run(ansatz=ansatz,qubit_hamiltonian=qubit_hamiltonian)