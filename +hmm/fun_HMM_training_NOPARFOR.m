%---------------
% HMM_training
%---------------
% 
% Computes t.p.m.,e.p.m.,LL calling hmmtrain and parameters in HmmParam
% 

function all_data=fun_HMM_training_NOPARFOR(sequence,gnunits,HmmParam)
%
% TOLERANCE
TOL=1e-8;
% 
% global screen
warning('OFF');
% parameters
% Mode=HmmParam.Mode;
VarStates=HmmParam.VarStates;
NumRuns=HmmParam.NumRuns;
NumSteps=HmmParam.NumSteps;
%
seq_tot=sequence;
ntrials=length(sequence);
numEmissions=gnunits+1;
% initialize variables for nested parfor loop
% linearize nested loops on VarStates and NumSteps
iterations=[NumSteps, numel(VarStates)];
NumIter=prod(iterations);
for ix=1:NumIter
    [j(ix),cnt(ix)]=ind2sub(iterations,ix);
end
temp_LL=zeros(NumIter,1); % collect data from the 10 steps
tmp=repmat(struct('tpm',[],'epm',[],'LL',[],'DProb',[],'trguess',[],'emisguess',[],...
    'logliks',[],'estemis',[],'esttr',[],'tol',[]),1,NumIter);
% parfor ix=1:NumIter
for ix=1:NumIter
    % parameters and initial guesses for emissions and transitions
    % run the training with 10 different values of the initial
    % TPM from 0.99 to 0.999 and pick the one with the highest
    % likelihood logliks
    % INITIALIZE TPM
    tmp(ix).DProb=0.990+0.01*(2*rand(VarStates(cnt(ix)),1)-1); % random diagonal entries in initial t.p.m. in [0.98,1]
    tmp(ix).trguess=zeros(VarStates(cnt(ix)));
    for ent=1:VarStates(cnt(ix))
        tmp(ix).trguess(ent,1:VarStates(cnt(ix)))=((1-tmp(ix).DProb(ent))/(VarStates(cnt(ix))-1))*ones(1,VarStates(cnt(ix)));          % initial guess for the estimated probability of transition from state i to state j(ix)
        tmp(ix).trguess(ent,ent)=tmp(ix).DProb(ent);      % the probability of transitions are 0.001
    end

    % INITIALIZE EPM
    % the rows of emisguess are the number of states. in each row, a column is
    % the probability that the i-th symbol is emitted. in our bernoulli case we
    % only have two symbols, no spike = 1 and spike = 2 (note we cannot put 0 and 1
    % instead because hmmtrain does not handle the symbol 0).
    tmp(ix).emisguess=rand(VarStates(cnt(ix)),gnunits+1);
    tmp(ix).tol=1;
    for it_cnt=1:NumRuns % run NumRuns times unless we reach tolerance tol<TOL
        % RANDOMIZE order of trial presentations for training to avoid
        % biases or getting stuck in a local min of
        % loglikel'd
        rand_trials=randperm(ntrials);
        % seq is a cell array of length cntrials, whose elements are the emission sequences
        % of each trial.
        seq=arrayfun(@(x)x.data,seq_tot(rand_trials(:)),'UniformOutput',false);
        if it_cnt==1 % first run uses initial guesses
            oldLL=0;
            oldGuessE=tmp(ix).emisguess;
            oldGuessTR=tmp(ix).trguess;               
            [tmp(ix).esttr,tmp(ix).estemis,tmp(ix).logliks]=hmm.hmmtrain(seq,tmp(ix).trguess,tmp(ix).emisguess,'Maxiterations',1);
        elseif it_cnt>1
            oldLL=tmp(ix).logliks;
            oldGuessE=tmp(ix).estemis;
            oldGuessTR=tmp(ix).esttr;
            [tmp(ix).esttr,tmp(ix).estemis,tmp(ix).logliks]=hmm.hmmtrain(seq,tmp(ix).esttr,tmp(ix).estemis,'Maxiterations',1);
        end
        % ITERATE BaumWelch until variation of LL, E,
        % TR is all less than tol or it_cnt>NumRuns
        tolLL=(abs(tmp(ix).logliks-oldLL)/(1+abs(oldLL)));
        tolTR=norm(tmp(ix).esttr - oldGuessTR,inf)/VarStates(cnt(ix));
        tolE=norm(tmp(ix).estemis - oldGuessE,inf)/numEmissions;
        tmp(ix).tol=max([tolLL tolE tolTR]);
        if tmp(ix).tol<TOL
            fprintf('Reached TOL=%0.03g->break\n',TOL);
            break
        end
    end
    tmp(ix).tpm=tmp(ix).esttr;
    tmp(ix).epm=tmp(ix).estemis;
    temp_LL(ix)=tmp(ix).logliks;
end
%-----------------------------------------------
% convert back to NumSteps and VarStates indices
all_data=repmat(struct('tpm',[],'epm',[],'LLtrain',[]),NumSteps,numel(VarStates));
for ix=1:NumIter
    all_data(j(ix),cnt(ix)).LLtrain=-2*temp_LL(ix); % collect data from the 10 steps
    all_data(j(ix),cnt(ix)).tpm=tmp(ix).tpm; % collect data from the 10 steps
    tmp(ix).epm(tmp(ix).epm==0)=eps; % regularize epm to avoid NaNs;
    all_data(j(ix),cnt(ix)).epm=tmp(ix).epm; % collect data from the 10 steps
end