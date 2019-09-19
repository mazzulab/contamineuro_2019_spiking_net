%% full HMM analysis

%----------------------------
% CREATE NEURAL ENSEMBLE
%----------------------------
% run this cell to create a new ensemble from the simulation
% 
if 0
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
    win_train=repmat(win,ntrials,1);
end

%%
% Instructions for fitting your own dataset
% 1) set 'if 0' on the first cell "create neural ensemble" on line 8
% 2) reformat your data into the 'spikes' structure array with dimension [ntrials,nunits] and field .spk containing spike times as columns array (in seconds)
% 3) create your own 'win_train' array with dimensions [ntrials, 2], where each row is the [start, end] times for each trial (spike times in 'spikes' must be consistently aligned with 'win_train')
%
%-------
% HMM
%-------
% model selection method
MODELSEL='XVAL'; % 'BIC';%'AIC';%
% HMM PARAMETERS
[ntrials, gnunits]=size(spikes);
hmmdir=fullfile('data','hmm'); if ~exist(hmmdir,'dir'); mkdir(hmmdir); end
filesave=fullfile(hmmdir,'HMM');
DATAIN=struct('spikes',spikes,'win',win_train,'METHOD',MODELSEL,'filesave',filesave); % experiment_animal_session same as before
res=hmm.funHMM(DATAIN);

% OUTPUT OF HMM FIT
% 'res' is a structure array with fields: 'hmm_data','hmm_bestfit','hmm_results',
%                                         'hmm_postfit','hmm_multispikes','HmmParam',
%                                         'win_train','colors','BestStateInd','LLtot'
% Important fields:
%     .hmm_bestfit.tpm: K x K transition probability matrix, where K is the number of hidden states
%     .hmm_bestfit.epm: K x (nunits+1) emission probability matrix, where K is the number of hidden states, the (n+1)-th column represents the probability of silence - you can safely drop it
%     .hmm_bestfit.LLtrain: -2*loglikelihood of the data 
% 
%     .hmm_results(i_trial).pStates: array of dim [K,time] with posterior probabilities of each state in trial i_trial
%     .hmm_results(i_trial).rates: array of dim [K,nunits] with local estimate of emissions (i.e., firing rates in each state) conditioned on observations in trial i_trial 
%     .hmm_results(i_trial).Logpseq: -2*loglikelihood from local observations in trial i_trial
%     
%     .hmm_postfit(i_trial).sequence: array of dimension [4,nseq] where columns represent detected states (intervals with prob(state)>0.8), in the order they appear in trial
%         i_trial, and rows represent state [onset,offset,duration,label].
%         
%      HmmParam.VarStates: number of states K chosen with model selection
%      HmmParam.BinSize: bin size (seconds)
%      HmmParam.MinDur: shortest duration of detected states; shorter states are excluded from hmm_postfit(i_trial).sequence
%      HmmParam.MinP: state detection probability 
%      HmmParam.NumSteps: number of independent EM runs used to non-convexity issues
%      HmmParam.NumRuns: maximum number of iterations each EM runs for
     
    
%% PLOTS
HmmParam=res.HmmParam;
HmmParam.VarStates=res.HmmParam.VarStates(res.BestStateInd);
hmm_bestfit=res.hmm_bestfit; hmm_results=res.hmm_results; hmm_postfit=res.hmm_postfit;
colors=aux.distinguishable_colors(max(HmmParam.VarStates,4));
% plot tpm and epm
hmm.plot_tpm_epm;
% plot trials
hmm.plot_trials;
