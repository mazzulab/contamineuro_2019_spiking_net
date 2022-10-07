function [LLtottemp,hmm_all_data,hmm_all_bestfit,temp_SkipSpikesSess]=fun_HMM_BIC_AIC(Spikes,win_train,HmmParam)

% BINNING SPIKES
% create sequences from 5 s prior to cue to 5 s post delivery
[sequence, temp_SkipSpikesSess]=hmm.fun_HMM_binning(Spikes,HmmParam,win_train);
[ntrials, gnunits]=size(Spikes);

LLtottemp=zeros(numel(HmmParam.VarStates),1);
T=0; % total number of datapoints
for trial=1:ntrials
    T=T+size(sequence(trial).data,2);
end
%---------
% TRAINING
%---------
hmm_all_data=hmm.fun_HMM_training(sequence,gnunits,HmmParam); 
% rows=NumSteps runs with different initial conditions
% cols=VarStates runs with different # of hidden states
% for each choice of number of states in VarStates, select the best out of
% NumSteps runs
% -> This is needed to avoid getting stuck in local minima (HMM is non-convex)
ind_step=zeros(1,numel(HmmParam.VarStates));
hmm_all_bestfit=repmat(struct('tpm',[],'epm',[],'LLtrain',[]),1,numel(HmmParam.VarStates));
for st_cnt=1:numel(HmmParam.VarStates)
    tempLL=cell2mat(arrayfun(@(x)x.LLtrain,hmm_all_data(:,st_cnt),'uniformoutput',false));
    [~,ind_step(st_cnt)]=min(tempLL); % find index of initial cond. with highest LL
    hmm_all_bestfit(1,st_cnt)=hmm_all_data(ind_step(st_cnt),st_cnt); % epm fit
end
LL=[hmm_all_bestfit(1,:).LLtrain]';
LLtottemp(1:numel(HmmParam.VarStates),1)=LL+HmmParam.NP(HmmParam.VarStates',gnunits,log(T));
