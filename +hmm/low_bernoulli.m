% This function turns a cell array of spike times for simultaneouly recorded
% units into a low-pass filtered matrix of Bernoulli emissions
% It eliminates the spikes occurring at negative times
% When more than one unit spike in the same bin, it randomly eliminates all
% but one (turns Poisson process into Bernoulli). 
% INPUT:  X{#units} = cell array, each element is the array of spike times of
%                    one unit
%         Tmax      = end of trial (s)
%         BinSize   = size of bin
%
% OUTPUT: matrix of emissions: 1st dim=bin; 2nd dim=unit
%         it is 0 for no spike, 1 for a spike

function [X_temp, varargout]=low_bernoulli(temp_shift,windiff,BinSize)
    
Bins=ceil(windiff/BinSize);
gnunits=numel(temp_shift);
X_temp=zeros(Bins,gnunits);
% eliminate pre-cue and post-delivery spikes
a_pos=cellfun(@(x)x(x(:)>0 & x(:)<=windiff),temp_shift,'UniformOutput',false);
% turn spike times into ms and ceil them to get indices of Bins
ind_spk=cellfun(@(x)ceil(x(:)/BinSize),a_pos,'UniformOutput',false);
% set X_temp indices to one for each spike
for unit=1:gnunits
    X_temp(ind_spk{unit}(:),unit)=1;
end
% if more than one unit spike in same bin, leave just one at
% random: turn poisson into bernoulli
%
% find bins with more than one spike
ind_many=find(sum(X_temp')>1);
% # of multiple spikes per bin in bins with more than one spike
SkipSpikes=sum(X_temp(ind_many,:),2);
for b=1:numel(ind_many)
    u=[];
    % find multiple units firing in that bin
    u=find(X_temp(ind_many(b),:));
    u_ind=randperm(numel(u));
    X_temp(ind_many(b),u(u_ind(1:end-1)))=0;
end
if nargout>0
    varargout={SkipSpikes};
end