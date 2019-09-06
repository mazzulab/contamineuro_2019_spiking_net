%% full HMM analysis

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
% model selection method
MODELSEL='XVAL'; % 'BIC';%'AIC';%
% HMM PARAMETERS
[ntrials, gnunits]=size(spikes);
win_train=repmat(win,ntrials,1);
hmmdir=fullfile('data','hmm'); if ~exist(hmmdir,'dir'); mkdir(hmmdir); end
filesave=fullfile(hmmdir,'HMM');
DATAIN=struct('spikes',spikes,'win',win_train,'METHOD',MODELSEL,'filesave',filesave); % experiment_animal_session same as before
res=hmm.funHMM(DATAIN);

%% PLOTS
HmmParam=res.HmmParam;
HmmParam.VarStates=res.HmmParam.VarStates(res.BestStateInd);
hmm_bestfit=res.hmm_bestfit; hmm_results=res.hmm_results; hmm_postfit=res.hmm_postfit;
colors=aux.distinguishable_colors(max(HmmParam.VarStates,4));
% plot tpm and epm
hmm.plot_tpm_epm;
% plot trials
hmm.plot_trials;
