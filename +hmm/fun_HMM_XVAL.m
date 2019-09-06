function [LLtotxval,hmm_all_data,temp_SkipSpikesSess]=fun_HMM_XVAL(Spikes,indxval,win_train,HmmParam)

%-------------------
% CROSS VALIDATION
%-------------------
% Leave L-out: training set has K=sum(ntrials)-L trials, test set has L trial.
% Perform K xval runs. 
% each HMM is run with 1 NumSteps only
%
K=numel(indxval);
[ntrials,gnunits]=size(Spikes);

[sequence, temp_SkipSpikesSess]=hmm.fun_HMM_binning(Spikes,HmmParam,win_train);
%
% LLtottemp=repmat(struct('z_out',[]),1,K); % z for each hold out and number of states
NStates=numel(HmmParam.VarStates);
LLtotxval=NaN(1,NStates); % store average LL of LLholdout values for each state

iterations=[K, numel(HmmParam.VarStates)];
NumIter=prod(iterations);
HmmTemp=repmat(HmmParam,1,NumIter);
jj=NaN(1,NumIter); cnt=NaN(1,NumIter);
for ix=1:NumIter
    [jj(ix),cnt(ix)]=ind2sub(iterations,ix);
    HmmTemp(ix).VarStates=HmmParam.VarStates(cnt(ix));
end
temp_LLholdout=NaN(1,NumIter);% store K holdout LL values for each state
temp_hmm_all_bestfit=repmat(struct('tpm',[],'epm',[],'LLtrain',[]),1,NumIter);
temp_results=repmat(struct('rates',[],'pStates',[],'Logpseq',[]),1,NumIter);
parfor ix=1:NumIter
        %---------
        % TRAINING
        %---------
        temp_hmm_all_bestfit(ix)=hmm.fun_HMM_training_NOPARFOR(sequence(indxval(jj(ix)).train),gnunits,HmmTemp(ix)); 
        % rows=NumSteps runs with different initial conditions
        % cols=VarStates runs with different # of hidden states
        %-----
        % HOLD-OUT
        %-----
        % HMM DECODING HOLD-OUT
        temp_results(ix).Logpseq=sum(arrayfun(@(x)[x(:).Logpseq],hmm.fun_HMM_decoding(Spikes(indxval(jj(ix)).test,:),...
            temp_hmm_all_bestfit(ix),HmmTemp(ix),win_train(indxval(jj(ix)).test,:))));
        temp_LLholdout(ix)=temp_results(ix).Logpseq; % LL of hold-out trials in k-th xval run
end
LLholdout=NaN(K,NStates);
hmm_all_data=repmat(struct('tpm',[],'epm',[],'LLtrain',[]),K,NStates);
for ix=1:NumIter
    LLholdout(jj(ix),cnt(ix))=temp_LLholdout(ix); % collect data from the 10 steps
    hmm_all_data(jj(ix),cnt(ix))=temp_hmm_all_bestfit(ix); % collect data from the 10 steps
end
LLtotxval(1,1:NStates)=mean(LLholdout,1);




