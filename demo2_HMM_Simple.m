%----------------------------
% CREATE NEURAL ENSEMBLE
%----------------------------
% select N neurons at random from the network to create a small ensemble for HMM fits
savedir=fullfile('data'); if ~exist(savedir,'dir'); mkdir(savedir); end % setup directory for saving HMM data
file_sim=fullfile(savedir,'results.mat');  % file where simulation results are saved
% option='random' -> select N neurons at random (eg, N=50)
% option='cluster' -> select N neurons per cluster (eg, N=2)
% option='random'; N=50;
option='cluster'; N=2;
[spikes,win]=aux.fun_create_ensemble(file_sim,option,N); 
% newspikes is a struct array with dimension [ntrials,nunits] and field .spk containing spike times
% windel=trial window

%%
%-------
% HMM
%-------
% HMM PARAMETERS
[ntrials, gnunits]=size(spikes);
win_train=repmat(win,ntrials,1);
HmmParam=struct();
% NUMBER OF HIDDEN STATES
HmmParam.VarStates=10; % number of hidden states
%--------------------
% HmmParam.AdjustT=0.1; % interval to skip at trial start to avoid canonical choice of 1st state in matlab
HmmParam.BinSize=0.002;%0.005; % time step of Markov chain
HmmParam.MinDur=0.05;   % min duration of an admissible state (s) in HMM DECODING
HmmParam.MinP=0.8;      % pstate>MinP for an admissible state in HMM ADMISSIBLE STATES
HmmParam.NumSteps=1;%    % 10 number of fits at fixed parameters to avoid non-convexity
HmmParam.NumRuns=500;%     % 50% % number of times we iterate hmmtrain over all trials

tic
% transform spike times into observation sequence
[sequence, ~]=hmm.fun_HMM_binning(spikes,HmmParam,win_train);
% train HMM
hmm_bestfit=hmm.fun_HMM_training_NOPARFOR(sequence,gnunits,HmmParam);
fprintf('HMM fit with %d states\n',HmmParam.VarStates);
toc
% estimate posterior probabilities
hmm_results=hmm.fun_HMM_decoding(spikes,hmm_bestfit,HmmParam,win_train);
% HMM ADMISSIBLE STATES -> state sequences
hmm_postfit=hmm.fun_HMM_postfit(spikes,hmm_results,HmmParam,win_train);

%% PLOTS
hmmdir=fullfile('data','hmm'); if ~exist(hmmdir,'dir'); mkdir(hmmdir); end
colors=aux.distinguishable_colors(max(HmmParam.VarStates,4));
% plot tpm and epm
hmm.plot_tpm_epm;
% plot trials
hmm.plot_trials;


