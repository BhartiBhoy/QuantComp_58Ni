---

## Shell-model study of <sup>58</sup>Ni  using quantum computing algorithm

authors:
  - name: Bharti Bhoy
    orcid: 0000-0003-0495-9916
    affiliation: 1 
  - name: Abhishek
    orcid: 0000-0003-2226-3146
    affiliation: 1
  - name: Paul Stevenson
    orcid: 0000-0003-2645-2569
    affiliation: 1
 
    
affiliations:
 - name: Department of Physics, University of Surrey, Guildford, Surrey, GU2 7XH, United Kingdom
   index: 1
bibliography: paper.bib

---


## Summary

This software package implements a quantum computing approach for solving the nuclear shell model problem, focusing on the <sup>58</sup>Ni  nucleus [@Bhoy_2024]. The software converts J-scheme two-body matrix elements into M-scheme matrix elements and maps nuclear orbitals onto qubits. It then uses these qubit-mapped elements to convert the Hamiltonian into the Jordan-Wigner (JW) transformation. By employing a Variational Quantum Eigensolver (VQE) with various classical optimizers, the software identifies the ground state energy and plots the corresponding wavefunction. The accuracy of the software is validated against classical shell model calculations, demonstrating its potential as a powerful tool for nuclear physics research.

## Statement of Need

In the realm of nuclear physics, traditional computational methods often struggle with the increasing complexity of simulating larger systems, such as mid-mass nuclei [@FGargano]. This challenge is particularly evident in the nuclear shell model, where the number of interacting particles significantly raises the computational demands, making exact solutions difficult to achieve with classical computers.

The software presented in this study addresses this critical need by leveraging quantum computing to perform shell model calculations on <sup>58</sup>Ni . By utilizing the Variational Quantum Eigensolver (VQE) and a problem-specific ansatz, this software provides an accurate simulation of the ground and low-lying excited states, reproducing exact energy values. This approach is scalable, versatile, and applicable to a broader range of nuclear systems, offering a new avenue for exploring nuclear structure that overcomes the limitations of classical methods.

Furthermore, this software serves as a benchmark for future quantum computing applications in nuclear physics, providing a vital comparison with classical shell model calculations. As quantum hardware advances, the methods implemented here will play a key role in extending the capabilities of nuclear simulations, making this software an essential tool for the research community.

## Acknowledgments

This work is supported by the UK STFC under the Quantum Technologies for Fundamental Physics programme, with Grant ST/W006472/1.

## References