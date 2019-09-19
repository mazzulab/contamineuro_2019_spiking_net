# contamineuro_2019_spiking_net
 Tutorial for Contamineuro Summer School 2019
 
Demo script for running spiking network simulations and analyses

by Luca Mazzucato 2019
----------------------------------------
Please cite:
L. Mazzucato, G. La Camera, A. Fontanini 
Expectation-induced modulation of metastable activity underlies faster coding of sensory stimuli, 
Nat. Neuro. 22, 787-796 (2019).
----------------------------------------

The tutorial demo1_simulation.m runs simulations of LIF clustered networks with excitatory (E) and inhibitory (I) spiking neurons.
You may run 2 different network architectures: 
1) A network with E clusters only (ClustersOption='E') [This part reproduces results from L. Mazzucato et al., 2019]
2) A network with E and I clusters (ClustersOption='EI')  [This part generates unpublished results (manuscript in preparation)]

The tutorial demo2_HMM_simple.m fits a Hidden Markov Model (HMM) with a fixed number of states to the network simulations from the previous scripts. This script is used to familiarize with HMM analyses.

The tutorial demo3_HMM_Full.m performs model selection for the number of HMM states, runs a full HMM fit and plots the results. It can be readily used on ensemble recordings from ephys data.

Scripts are optimized for parallel computation on multi-core machines.