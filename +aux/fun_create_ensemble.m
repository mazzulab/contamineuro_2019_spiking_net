% function LIF2_PreProcess(SessSet,name,Opt)
% RUN BEFORE LIF2_HMM
% Load simulation files -> load firings var -> creates spike times structure
% Call ElimSilent to eliminate units <2 spk/s; create ensembles; save

function [newspikes, win]=fun_create_ensemble(file_sim,option,N)

hmmdir=fullfile('data'); if ~exist(hmmdir,'dir'); mkdir(hmmdir); end % setup directory for saving HMM data
file_ensemble=fullfile(hmmdir,'ensemble.mat'); % file were ensemble spikes are saved

dataload=load(file_sim);
params=dataload.params;
firings=dataload.firings;
aux.v2struct(params);

% CREATE SPIKE STRUCTURE
% 1) pick NUnits units from 'firings' structure
% 2) eliminate silent units
tic
% load(paramsfile,'Sim','N_e','p');
TotUnits=f*p*N_e;% # of clustered neurons
win=[Sim.t_Start,Sim.t_End];
ind_units=1:TotUnits;
% extract NUnits random exc units from 'firings'
NUnits=min(N_e,5*max(N)); % 5x the largest ensemble
temp_firing=firings;
ntrials=numel(temp_firing);
spikes=[]; spikes(ntrials,NUnits).spk=[];
for trial=1:ntrials
    for unit=1:numel(ind_units)
        ind=(temp_firing{trial}(:,2)==ind_units(unit));
        spikes(trial,unit).spk=temp_firing{trial}(ind,1);
    end
end

if strcmp(option,'random')
    % Eliminate silent units
    spikes=aux.ElimSilent(spikes,win);
    % Recheck
    spikes=aux.RecheckSilent(spikes,win);
    a=randperm(size(spikes,2));
    temp_spikes=spikes(:,a(1:N)); % all newspikes
    fprintf('\n            picked %d units\n',N);
elseif strcmp(option,'cluster')
    
    
    NcE=params.popsize';
    % find cluster indices from weight matrix or file
    MinCluster=N;% largest # of neurons per cluster needed
    % Extract spikes  & find silent clusters
    indStart=[0 cumsum(NcE)];
    ind_units=repmat(struct('ind',[]),1,p);
    for cl=1:p
        Results(cl).tempSpikes=[];
        % without overlaps
        ind_units(cl).ind=indStart(cl)+1:indStart(cl+1);
        % with overlaps
    end
    for cl=1:p
        temp_firing=firings;
        ntrials=numel(temp_firing);
        spikes=[]; spikes(ntrials,numel(ind_units(cl).ind)).spk=[];
        for trial=1:ntrials
            for unit=1:numel(ind_units(cl).ind)
                ind=(temp_firing{trial}(:,2)==ind_units(cl).ind(unit));
                spikes(trial,unit).spk=temp_firing{trial}(ind,1);
            end
        end
        Results(cl).tempSpikes=spikes;
        Results(cl).tempSpikes=aux.ElimSilent(Results(cl).tempSpikes,win);
    end
    ClusterIndex=[];
    for cl=1:p
        if isempty(Results(cl).tempSpikes)
            ClusterIndex=[ClusterIndex cl]; % remove silent cluster
            fprintf('Cluster %d inactive, skip it...\n',cl);
        elseif size(Results(cl).tempSpikes,2)<MinCluster
            ClusterIndex=[ClusterIndex cl]; % remove silent cluster
        end
    end
    ClusterIndex=setxor(ClusterIndex,1:p);
    fprintf('%d active clusters\n',numel(ClusterIndex));
    temp_spikes=[];
    for cl=1:numel(ClusterIndex)
        temp_newspikes=Results(cl).tempSpikes;%sim_spikes(cl);
        fprintf('     cluster #%d:\n',cl);
        if size(temp_newspikes,2)<N
            fprintf('     ... Failed to find %d non silent unit in cluster %d\n',N,cl);
        else
            atemp=1:size(temp_newspikes,2);
            transfer_spikes=temp_newspikes(:,atemp(1:N));
            temp_spikes=[temp_spikes transfer_spikes]; % all newspikes
            fprintf('            picked %d units\n',N);
        end
    end
end
newspikes=temp_spikes;

save(file_ensemble,'newspikes','params');
