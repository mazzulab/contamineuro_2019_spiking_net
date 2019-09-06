function temp_newspikes=RecheckSilent(temp_newspikes,windel)

ntrials=size(temp_newspikes,1);
%
ind=[];
bins=[windel(1) 0 windel(2)];
tempRateBins=Spikes2Bins(temp_newspikes,bins);
for b=numel(bins)-1
    tempRate=squeeze(tempRateBins(:,b,:));
    lambda=mean(tempRate);
    % collect units with low rates
    ind(b).Zero=find(lambda<1);
end
% eliminate neurons that are silent during all evoked events or all spontaneous
% events
indZero=[];
for b=1:numel(bins)-1
    indZero=[indZero ind(b).Zero];
end
if ~isempty(indZero)
    indRemove=unique(indZero);
    temp_newspikes=temp_newspikes...
        (1:ntrials,setxor(indRemove,1:size(temp_newspikes,2)));
end