% function Spikes2Seq(spikes,windel,BinSize)
% 
% Convert ensemble spike train spikes(:).spk in the interval windel 
% into a sequence of observations seq in bins of BinSize s
% 
% Luca Mazzucato 2015

function [seq, varargout]=Spikes2Seq(spikes,windel,BinSize)
        
% create cell array with spike times for all units in this
% trial
% create cell array with spike times for all units in this
% trial, shifting by win(1)
temp_all=arrayfun(@(x)x.spk,spikes,'UniformOutput',false);
temp_shift=cellfun(@(x)x-windel(1),temp_all(:),'UniformOutput',false);
[X_temp, SkipSpikes]=hmm.low_bernoulli(temp_shift',windel(2)-windel(1),BinSize);
if nargout>0
    varargout={SkipSpikes};
end
%------------------
% TURN INTO HMM FORMAT
%------------------
Bins=round((windel(2)-windel(1))/BinSize);
%         fprintf('\n%d Bins \n',Bins);
Xseq=X_temp'; % rows=units; cols=bins;
[i, j]=ind2sub(size(Xseq),find(Xseq>0));
seq=(size(Xseq,1)+1)*ones(1,Bins);
seq(j)=i;
% format ready for hmmtrain
% length=Bins+1; each entry is an emission with following
% values
% gnunits+1=no spikes (default)
% i=spike fired by unit i
% 
%             fprintf('\n length of seq=%d\n',numel(seq));
