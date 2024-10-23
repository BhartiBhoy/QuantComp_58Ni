# Shell-model study of $^{58}\text{Ni}$ using quantum computing algorithm

This repository contains the code used for the study of energy levels in the 58Ni nucleus, using quantum computing techniques. Specifically, we utilize the Variational Quantum Eigensolver (VQE) algorithm in combination with a problem-specific ansatz. This study demonstrates how quantum algorithms can be applied to nuclear structure studies and compares results with classical shell model computations.

## Objective
The primary goal is to calculate the low-lying energy levels (ground and excited states) of 58Ni, a nucleus relevant to nuclear physics and astrophysical reactions. We aim to demonstrate that quantum computing techniques can achieve results comparable to classical methods, thus pushing the boundaries of quantum applications in nuclear physics.

## Methodology
- **Quantum Algorithm**: We use the VQE, a hybrid quantum-classical algorithm, to minimize the energy expectation value of the 58Ni Hamiltonian and compute its ground and excited state energies.
- **Ansatz**: The wavefunction ansatz is specifically tailored for the 58Ni nucleus to represent the desired quantum states (ground, first, and second excited states).
- **Qubit Mapping**: Fermionic creation and annihilation operators are mapped onto qubits using the Jordan-Wigner transformation. Each qubit corresponds to a particular nuclear state in the shell model.
- **Optimization**: Optimization algorithms like COBYLA, SLSQP, SPSA and Gradient-Descent are used to adjust the ansatz parameters and minimize the Hamiltonian.


## Results
- The VQE simulation accurately reproduces the ground state and first and second excited state energies of 58Ni.
- Quantum results match well with classical shell model diagonalization, confirming the validity of the quantum approach.
- The study demonstrates the potential of quantum computing in nuclear physics, particularly in nuclear structure studies.

## Dependencies
- Python 3.11
- Qiskit
- Numpy
- Matplotlib
- Scienceplots

## How to Run
1. Clone the repository:
   ```bash
   git clone https://github.com/BhartiBhoy/58Ni_SM_VQE.git
   ```
2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the main.py file in the root directory:
   ```bash
   python main.py
   ```
   By default it will calculate for the ground state and will use the COBYLA optimizer. 
   To select a different excited state and different optimizer (COBYLA, SPSA, SLSQP, GradientDescent), there are flags that can be used as follows:
   ```bash
   python main.py --state 1 --optimizer GradientDescent
   or
   python main.py -s 1 -o GradientDescent
   ```
   This will select the first excited state (-s 1) and will user the GradientDescent optimizer. 
   The second excited state can be selected using -s 2.


## Citation
If you use this code in your research, please cite the following paper:
```
Bharti Bhoy and Paul Stevenson, "Shell-model study of 58Ni using quantum computing algorithm", New J. Phys. 26 (2024), 075001. DOI: 10.1088/1367-2630/ad5756
```

## License
This project is licensed under the Creative Commons Attribution 4.0 International License. 
