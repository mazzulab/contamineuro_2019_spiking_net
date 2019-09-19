
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
% spikes is a struct array with dimension [ntrials,nunits] and field .spk containing spike times
% windel=trial window

%%
% Instructions for fitting your own dataset
% 1) Comment the first cell "create neural ensemble"
% 2) reformat your data into the 'spikes' structure array with dimension [ntrials,nunits] and field .spk containing spike times as columns array (in seconds)
% 3) create your own 'win_train' array with dimensions [ntrials, 2], where each row is the [start, end] times for each trial (spike times in 'spikes' must be consistently aligned with 'win_train')
%-------
% HMM
%-------
% HMM PARAMETERS
[ntrials, gnunits]=size(spikes);
win_train=repmat(win,ntrials,1);
HmmParam=struct();
% NUMBER OF HIDDEN STATES
HmmParam.VarStates=5; % number of hidden states
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

% OUTPUT OF HMM FIT
%
% Important variables:
%     hmm_bestfit.tpm: K x K transition probability matrix, where K is the number of hidden states
%     hmm_bestfit.epm: K x (nunits+1) emission probability matrix, where K is the number of hidden states, the (n+1)-th column represents the probability of silence - you can safely drop it
%     hmm_bestfit.LLtrain: -2*loglikelihood of the data 
% 
%     hmm_results(i_trial).pStates: array of dim [K,time] with posterior probabilities of each state in trial i_trial
%     hmm_results(i_trial).rates: array of dim [K,nunits] with local estimate of emissions (i.e., firing rates in each state) conditioned on observations in trial i_trial 
%     hmm_results(i_trial).Logpseq: -2*loglikelihood from local observations in trial i_trial
%     
%     hmm_postfit(i_trial).sequence: array of dimension [4,nseq] where columns represent detected states (intervals with prob(state)>0.8), in the order they appear in trial
%         i_trial, and rows represent state [onset,offset,duration,label].


%% PLOTS
hmmdir=fullfile('data','hmm'); if ~exist(hmmdir,'dir'); mkdir(hmmdir); end
colors=aux.distinguishable_colors(max(HmmParam.VarStates,4));
% plot tpm and epm
hmm.plot_tpm_epm;
% plot trials
hmm.plot_trials;


