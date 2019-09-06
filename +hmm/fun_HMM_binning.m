%----------------
% HMM_binning
%----------------
%
% varargout{1}=TotSkipSpikes.(trig{E})(trial).mspikes, array containing
% number of multiple spikes per bin, in bins with more than one spike per
% bin

function [sequence, TotSkipSpikes]=fun_HMM_binning(spikes,HmmParam,win_train,varargin)

if ~isempty(varargin)
    funBin=str2func(varargin{1});
else
    funBin=@hmm.Spikes2Seq;
end
warning('off');
BinSize=HmmParam.BinSize;
[ntrials]=size(spikes,1);
sequence=repmat(struct('data',[]),1,ntrials);
TotSkipSpikes=repmat(struct('mspikes',[]),1,ntrials);
for trial=1:ntrials
    win=win_train(trial,:);
    [seq, SkipSpikes]=funBin(spikes(trial,:),win,BinSize);
    sequence(trial).data=seq;
    TotSkipSpikes(trial).mspikes=SkipSpikes;
end
