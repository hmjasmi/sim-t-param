# T-Parameters Based Modeling for Stacked Intelligent Metasurfaces: Tractable and Physically Consistent Model

This repository contains the code and data associated with the paper:

H. Yahya, M. Nerini, B. Clerckx and M. Debbah, "[T-Parameters Based Modeling for Stacked Intelligent Metasurfaces: Tractable and Physically Consistent Model](https://arxiv.org/)," Archive, 2025.

## Table of Contents
- [Introduction](#introduction)
- [Code Overview](#code-overview)
- [Routines](#routines)
- [Results](#results)
- [Citation](#citation)

## Introduction
This repository provides the source code for the T-prameters based design presented in the paper. The goal is to demonstrate the effectiveness of the T-prameters based modeling to optimize SIM.

## Code Overview
The provided code includes the following MATLAB scripts:
- **func_MRT_GC.m**: A function that designs RIS scattering matrix to achieve ([passive MRT](https://ieeexplore.ieee.org/abstract/document/10771739)).
- **func_S2T.m**: A function that converts an S-parameters matrix to an T-parameters matrix.
- **func_MRT_init.m**: A function that provides the SIM phase shifts initialization based on the passive MRT and simplified channel.
- **func_sR_MAX_GDA_ExactExact.m**: A function that optimizes the SIM phase shifts to maximize the sum rate using GDA. The optimization and evaluation are based on exact channel.
- **func_sR_MAX_GDA_SimplExact.m**: A function that optimizes the SIM phase shifts to maximize the sum rate using GDA. The optimization is based on simplified channel and evaluation is based on {exact, simplified} channel.
- **func_SIM_MC_dipole.m**: A function that computes \mathbf{Z}_l for SIM based on ([dipole mutual impedances expression](https://ieeexplore.ieee.org/abstract/document/9319694)).  
- **func_SIM_RaySom.m**: A function that computes \mathbf{S}_{l,21} for SIM based on ([Rayleigh Sommerfeld diffraction coefficients](https://ieeexplore.ieee.org/abstract/document/10279173)). 

### Dependencies
- MATLAB

## Routines
The routine files are used to generate the results. The raw results are stored in Results folder.

## Results
The results are plotted and stored in Fig_saved folder.

## Citation
If you find this work useful in your research, please consider citing our paper:

@article{yahya2025sim-t-param,
  title={T-Parameters Based Modeling for Stacked Intelligent Metasurfaces: Tractable and Physically Consistent Model},
  author={Yahya, Hamad and Nerini, Matteo and Clerckx, Bruno and Debbah, Merouane},
  journal={arXiv preprint arXiv:},
  year={2025}
}
