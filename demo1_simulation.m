% Demo script for running spiking network simulations and analyses
% 
% by Luca Mazzucato 2019
%
% ----------------------------------------
% Please cite:
% L. Mazzucato, G. La Camera, A. Fontanini 
% Expectation-induced modulation of metastable activity underlies faster coding of sensory stimuli, 
% Nat. Neuro. 22, 787-796 (2019).
% ----------------------------------------
% This script run simulations of LIF networks with excitatory (E) and inhibitory (I) spiking neurons.
% You may run 2 different network architectures: 
% 1) A network with E clusters only (ClustersOption='E') [This part reproduces results from L. Mazzucato et al., 2019]
% 2) A network with E and I clusters (ClustersOption='EI')  [This part generates unpublished results (manuscript in preparation)]
%
%% CHOOSE CLUSTERED ARCHITECTURE
%------------------------
% ClustersOption='EI';%
ClustersOption='E';%
%------------------------
% LOAD PARAMETERS
%------------------------
paramsfile='params.mat'; % file where all network parameters are saved
if strcmp(ClustersOption,'E'); aux.create_params(paramsfile);
elseif strcmp(ClustersOption,'EI'); aux.create_params_EI(paramsfile);
end

%% CHOOSE STIMULI (OR CREATE YOUR OWN)
% You have 3 built-in stimulus templates:
% 0) no stimuli: run the network during ongoing activity
% 1) 'US' stimulus onset at t=0s (stimulus is on until the end of the trial);
%   a linearly ramping external current delivered to 50% of E clusters (chosen randomly, 
%   half of the neurons in each stimulus-selective clusters receive stimulus) with a slope of 
%   gain=0.2*mu_ext, where mu_ext is the baseline external current (bias)
% 2) 'CSgauss' stimulus onset at t=-0.5s (stimulus is on until the end of
%   the trial); double exponential profile with rise and decay times
%   [0.5,1]s; the stimulus targets all E neurons. For each neuron, the
%   stimulus peak is drawn from a gaussian distribution with mean 0 and
%   standard deviation 0.2*mu_ext.
% CUSTOM OPTIONS:
% 1) edit your custom-made stimuli inside aux.create_params. Available options: 'US', 'CSgauss' from L. Mazzucato et al., Nat. Neuro. 2019
% 2) list which stimuli to run in current trials in the cell array 'stimuli'
% 3) if stimuli={}, run trial with spontaneous activity only, without any stimulus
%------------------------
% NOTE:
% This demo runs the unexpected (UT) and expected (ET) conditions from L. Mazzucato et al., 2019, depending on the following settings:
% % 1) with condition='UT' you simulate an unexpected trial, with a taste
% stimulus delivered at t=0s and not anticipatory cue (Fig. 3a right panel in the paper)
% 2) with condition='ET' you simulate an expected trial, with a taste stimulus
% delivered at t=0, preceded by an anticipatory cue at t=-0.5s (Fig. 3a left panel in the paper).
%-----------------
% SELECT STIMULUS
%-----------------
% stimuli={}; % ongoing activity
% stimuli={'US'}; % stimulus evoked-activity targeting selective clusters
% stimuli={'CSgauss'}; % anticipatory cue speeds up network dynamics
stimuli={'US','CSgauss'}; % anticipatory cue preceeds stimulu delivery
save(paramsfile,'stimuli','-append');
savedir=fullfile('data'); if ~exist(savedir,'dir'); mkdir(savedir); end % setup directory for saving HMM data

%% RUN SIMULATION
ntrials=20; % number of trials
file_sim=fullfile(savedir,'results.mat');  % file where simulation results are saved
%---------------------------
% GENERATE SYNAPTIC WEIGHTS
%---------------------------
% J = N x N matrix of synaptic weights
% params = structure containing all network parameters
if strcmp(ClustersOption,'E'); [J, params]=aux.fun_SynWeights(paramsfile);
elseif strcmp(ClustersOption,'EI'); [J, params]=aux.fun_SynWeights_EI(paramsfile);
end
[stimulus_save, params]=aux.fun_stim(params); % STIMULUS
%------------------------
% SIMULATION
%------------------------
tic
firings=cell(1,ntrials); % cell array with all spike times in each trial
PlotData=cell(1,ntrials); % cell array with data for plotting
% parfor iTrial=1:ntrials % uncomment this line if you have a multi-core  machine with 4 or more cores
for iTrial=1:ntrials
    ParamsRun=params;
    ParamsRun.Ext=stimulus_save.Ext;
    ParamsRun.Stimulus=stimulus_save.Stimulus;
    ParamsRun.J=J;
    fprintf('--- Start SIM ...\n');
    [firings{iTrial}, PlotData{iTrial}]=aux.fun_LIF_SIM(ParamsRun);
end
% SAVE results
save(file_sim,'params','firings','PlotData','stimulus_save');
fprintf('\nDone. Simulation saved in %s\n',file_sim);
toc
%%
%------------------------
% PLOT EVENTS
%------------------------
iTrial=1; % pick which trial to plot
dataload=load(file_sim); % load simulation results
data=dataload.PlotData{iTrial}; % membrane potentials
firings=dataload.firings{iTrial}; % spikes
Params=dataload.params; % parameters
Params.Ext=dataload.stimulus_save.Ext; % external currents
Params.savedir=savedir;
aux.fun_PlotTrial(data,firings,Params);
% figure 1 - rasterplot of all neurons in trial
% figure 2 - time course of membrane potential and PSC traces for E and I representative neurons
% figure 3 - time course of firing rate in clusters
% figure 4 - time course of stimuli (with CSgauss stimulus, the cue profile should be multiplied by a factor drawn from figure 5, one for each neuron
% figure 5 (with CSgauss stimulus only) - across-neurons distribution of cue peak values
