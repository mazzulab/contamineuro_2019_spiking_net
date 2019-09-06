%---------------------
% HMM_postfit
%---------------------

% POST FIT PROCESSING

function hmm_all_postfit=fun_HMM_postfit(Spikes,hmm_all_results,HmmParam,win_train)
 
% ADMISSIBLE states = states for which p>HmmParam.MinP for at least MinDur consecutive bins
%
% 1) keep only states whose probability pstates>HmmParam.MinP for at least MinDur ms ("admissible states"):
% 2) each admissible state is a point in the gnunits-dimensional
% space on which we run Kmeans
% 
% hmm_postfit(ntrials).sequence  
% 1st row: start of admissible states in sequence
% 2nd row: end of admissible states in sequence
% 3rd row: duration (end-start)
% 4th row: state number for each event (same as in t.p.m.)
[ntrials,~]=size(Spikes);
hmm_all_postfit=repmat(struct('sequence',[]),1,ntrials);
for trial=1:ntrials
    win_coding=win_train(trial,:);
    pstates_rate=hmm_all_results(trial).pStates;
%     numstates=size(pstates_rate,1);
    temp_sequence=hmm.fun_HMM_postfit_sequence(pstates_rate,HmmParam,win_coding); % script containing the state segmentation
    hmm_all_postfit(trial).sequence=temp_sequence;
end
