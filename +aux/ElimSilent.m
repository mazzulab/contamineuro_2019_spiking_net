%----------------------------------------------------------------

% function newspikes=ElimSilent(sim_spikes,windel,varargin)
%
% Eliminate units firing at less than 1 Hz and selects number of units for
% creating ensembles (may load same numbers as in data)
% INPUT: structure sim_spikes(NTrials,NUnits).spk=spikes times
%        windel=[t_start t_end]: trial window (s)
%        varargin{1}=if numeric, it is the number of output units in the ensemble
% OUTPUT: structure newspikes(NTrials,MaxUnits).spk=spikes times
%                 where MaxUnits=varargin{1} or original size


function newspikes=ElimSilent(sim_spikes,windel,varargin)

BinSilent=0.2;
temp=sim_spikes;
SimUnits=size(temp,2);
%
% ENSEMBLE OPTIONS
if ~isempty(varargin)
    method=varargin{1};
    if isnumeric(method)
        if method>=SimUnits
            fprintf('--- SIM ensembles have %d units, keeping that...\n',SimUnits);
            MaxUnits=SimUnits;
        else
            MaxUnits=method;
            fprintf('--- ElimSilent: Picking %d units/ensemble...\n',method);
        end
    else
        fprintf('\n--- Error in optional argument, see help...\n');
    end
else
    MaxUnits=SimUnits;
end


newspikes=[];

spikecnt=[];
nunits=size(sim_spikes,2);
spikecnt(nunits).spk=[];
% SPIKECOUNT
spks=sim_spikes;
[ntrials, nunits]=size(spks);
bins=windel(1):BinSilent:windel(2);
nbins=length(bins)-1; % number of bins
for unit=1:nunits
    spcnt_bis=[];
    for b=2:length(bins)
        spcnt=[];
        for trial=1:ntrials
            sp=numel(find(spks(trial,unit).spk<bins(b) & spks(trial,unit).spk>=bins(b-1))); % spikes per bin
            spcnt=[spcnt;sp];
        end
        spcnt_bis=[spcnt_bis spcnt];
    end
    spikecnt(unit).spk=[spikecnt(unit).spk; spcnt_bis]; % row are trials, columns are bins
end
% GOOD UNITS
% eliminate neurons whose average rate is less than 2 Hz in all
% events.
% gu is the set of good neurons. run units loops using gu.
gu=[];
% FIX # OF UNITS PER ENSEMBLE
u_cnt=0;
UnitSet=randperm(nunits);
for unit=1:numel(UnitSet)
    if u_cnt<MaxUnits
        if  mean(mean(spikecnt(UnitSet(unit)).spk,1))/BinSilent>1
            gu=[gu,UnitSet(unit)];
            u_cnt=u_cnt+1;
        end
    else
        break;
    end
end
gu=unique(gu);
% REPLACE spikes and spikecount with newspikes and gucount BY KEEPING ONLY GOOD UNITS GU
% skip if there are no good units
if any(gu)
    for unit=1:numel(gu)
        for trial=1:ntrials
            newspikes(trial,unit).spk=sim_spikes(trial,gu(unit)).spk;
        end
    end
end
