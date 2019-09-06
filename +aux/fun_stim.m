function [stimulus_save, params]=fun_stim(params)

% unpack vars
p=params.p;
Stimulus=params.Stimulus;
events=params.events;
Sim=params.Sim;

% for each event, create external current and  stim properties in .Ext 
% and add clusters selective to the event .Stimulus.feat(n).StimClust
stimulus_save=struct('Ext',[],'Stimulus',[]);

% select stimuli
stimuli=events;
temp_Stimulus=struct('input',Stimulus.input);
indfeat=zeros(1,numel(stimuli));
for ev=1:numel(stimuli)
    fprintf('Stimulus %s',stimuli{ev});
    % match current stimuli to features in Stimulus
    indfeat(ev)=find(cell2mat(arrayfun(@(x)strcmp(stimuli{ev},x(:).name),...
        Stimulus.feat(:),'uniformoutput',false)));
end
fprintf('\n');
if ~isempty(indfeat)
    temp_Stimulus.feat(1:numel(indfeat))=Stimulus.feat(indfeat);
    for n=1:numel(indfeat)
        sclust=[];
        if ~isempty(temp_Stimulus.feat(n).selective)
            sclust=find(temp_Stimulus.feat(n).selective(1,:));
        end
        temp_Stimulus.feat(n).StimClust=sclust;
    end
end
Stimulus=temp_Stimulus;
Ext=struct('Mu',[]);

% LOAD PARAMETERS
fieldNames={'Sim','Network','p','popsize','clustermatrix','N_e','N_i','Cext','Jee_ext','Jie_ext','ni_ext','tau_e','tau_i','fieldNames'};
aux.v2struct(params,fieldNames);
cusumNcE=[0 cumsum(popsize)'];
Tseq=Sim.t_Start:Sim.dt_step:Sim.t_End; 

if ~isempty(events)
    feat=Stimulus.feat;
    nstim=numel(feat); % number of stimuli in current trials
    stim=repmat(struct('profile',[],'ind',[],'interval',[]),1,nstim);
    temp_ind=repmat(struct('ind',[]),1,nstim); % stores indices for mixed cue (see below)
    for n=1:nstim
        % stimulus interval
        interv=feat(n).interval;
        Gain=feat(n).gain;
        if ~isempty(strfind(feat(n).name,'gauss'))
            Gain=1; % with gaussian stim set profile to peak at 1, then multiply each profile by gaussian with SD feat(n).gain for each neuron in feat(n).gauss
        end
        Profile=feat(n).profile;
        Profile=@(t)Profile(t-interv(1));
        MaxThInput=max(abs(Profile(Tseq(Tseq>interv(1) & Tseq<interv(2)))));
        Profile=@(t)Gain*Profile(t)/MaxThInput;
        stim(n).profile=@(t)Profile(t); % fraction increase above baseline
        % selective neurons
        StimClust=Stimulus.feat(n).StimClust; % clusters activated by current stimulus
        % units selective to stimulus
        ind=[]; % indices of stim sel units
        switch feat(n).selectivity
            case 'mixed'
                for c=StimClust
                    pop_ind=find(clustermatrix(:,c));
                    for k=1:numel(pop_ind)
                        ind=[ind cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                    end
                end
            case 'exc'
                ind=1:N_e;
            otherwise
                ind=1:N_e;
        end
        % sparsify
        a=randperm(numel(ind));
        temp_ind(n).ind=ind;
        ind=ind(a(1:round(feat(n).connectivity*numel(ind))));
        % gaussian stimulus, draw from randn
        if ~isempty(strfind(feat(n).name,'gauss'))
            stim(n).gauss=feat(n).gain*randn(numel(ind),1);
        end
        %
        stim(n).ind=ind;
        stim(n).interval=interv;
        stim(n).name=feat(n).name;
        stim(n).StimClust=StimClust;
        stim(n).selectivity=feat(n).selectivity;
        
    end
    Ext.stim=stim;
end 
Ext.Mu=params.Mu;



stimulus_save.Ext=Ext;
stimulus_save.Stimulus=temp_Stimulus;



