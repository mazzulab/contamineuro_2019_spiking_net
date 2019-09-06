% SIM of one trials given parameters
%
% Luca Mazzucato March 2014

% SET OPTIONS
% ParamsRun = structure containing parameters for simulation


function [all_firings, PlotData]=fun_LIF_SIM(ParamsRun)


Theta=ParamsRun.Theta; Sim=ParamsRun.Sim; events=ParamsRun.events;%Stimulus=ParamsRun.Stimulus;
Ext=ParamsRun.Ext; J=ParamsRun.J; N_e=ParamsRun.N_e; N_i=ParamsRun.N_i; p=ParamsRun.p; He=ParamsRun.He;
Hi=ParamsRun.Hi; tau_e=ParamsRun.tau_e; tau_i=ParamsRun.tau_i; tausyn_e=ParamsRun.tausyn_e;
tausyn_i=ParamsRun.tausyn_i; tau_arp=ParamsRun.tau_arp;
%
all_firings=[];
dt=Sim.dt_step;            % time step (s)
Tseq=Sim.t_Start:dt:Sim.t_End;  

%--------------------
% PARAMETERS
%--------------------
% CELL
VEreset=He*Theta(1);
VIreset=Hi*Theta(end);
%
%----------
% STIMULUS
%----------
% build external current
% baseline current for all neurons
% add stimuli on top of baseline: for each stimulus provide 
%              - profile (perc. increase on baseline current)
%              - index of selective neurons
%
% BASELINE EXTERNAL CURRENT
mu=Ext.Mu; % mu is an (N_e+N_i)-th array
if ~isempty(events)
    stim=Ext.stim;
    nstim=numel(stim); % number of stimuli in current trials
end
%----------------
% SYNAPTIC FILTER
%----------------

Tau.tausyn_e=tausyn_e; % exc synaptic time (fall time)
Tau.tausyn_i=tausyn_i; % exc synaptic time (fall time)
F=synaptic_trace(Tau,dt,N_e,N_i); % traces for recurrent connections

%--------------------------
% SIMULATION
%--------------------------
% preallocate memory for stored variable firings_tmp
% INITIAL CONDITIONS: random
v=[(Theta(1)-VEreset)/2*ones(N_e,1)+(Theta(1)-VEreset)/2*(2*rand(N_e,1)-1);...
    (Theta(end)-VIreset)/2*ones(N_i,1)+(Theta(end)-VIreset)/2*(2*rand(N_i,1)-1)];  % Initial values of v
% THRESHOLD VECTOR
VTh=[Theta(1)*ones(N_e,1); Theta(end)*ones(N_i,1)];
c=[VEreset*ones(N_e,1);  VIreset*ones(N_i,1)]; % reset potentials
% fprintf('\nVEThres=%g --- VIThres=%g',VTh(1),VTh(end));
% fprintf('\nVEreset=%g --- VIreset=%g \n',c(1),c(end));
% Excitatory neurons        Inhibitory neurons
tau=[tau_e*ones(N_e,1);       tau_i*ones(N_i,1)];
%
firings=zeros(10*numel(Tseq),2);
firings_cnt=0;
tic
%--------------------
% PLOT
%--------------------
PlotData=[];
PlotData.Ne_plot=N_e; % number of exc neuron to plot
PlotData.Ni_plot=N_i; % number of inh neurons to plot
ind_plot=[5; N_e+5]; % indices of neurons to plot
if ~isempty(events)
    indcue=find(cellfun(@(x)~isempty(x),strfind({stim(:).name},'CS')));
    if ~isempty(indcue)
        ind_plot(1)=stim(indcue).ind(1);
    end
end
nplot=numel(ind_plot); % number of neurons to plot (membrane potential plot)
vi=0; % running index for vplot
PlotData.vplot = zeros(nplot,round(Sim.plot_length/dt)); % store membrane potential for plots; rows=neurons, cols=time steps;
PlotData.iEplot = zeros(2,round(Sim.plot_length/dt)); % store EPSC for plots; rows=neurons, cols=time steps;
PlotData.iExtplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.iIplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.p=p;
PlotData.VTh=VTh;
PlotData.tau=tau;
PlotData.ind_plot=ind_plot;
%----------------------------
% RUN
%----------------------------

refr=zeros(size(mu,1),1);       % neurons in refractory state
for t=1:numel(Tseq)         % siMulation of 1000 ms 
    fired=find(v>VTh); % indices of spikes
    Isyn=zeros(N_e+N_i,1);
    % spikes
    if ~isempty(fired)      
        v(fired)=c(fired);  
        refr(fired)=tau_arp;
    end
    % recurrent synaptic current
    F=syn_evolve(F,fired);    
    % integrate
    muRun=mu;
    if ~isempty(events)
        for n=1:nstim
            if Tseq(t)>=stim(n).interval(1) && Tseq(t)<=stim(n).interval(2) 
                if strcmp(stim(n).name,'CSgauss')
                    muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
                else
                    muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind);
                end
            end
        end
    end
    Isyn=Isyn+J*F.f;
    v=v-v*dt./tau+muRun(:)*dt+Isyn*dt;
    % neurons in refractory state
    refr=max(-0.001,refr-dt); 
    v(refr>0)=c(refr>0);
    % store spikes
    if ~isempty(fired)
        % if firings_tmp has no more space, preallocate more memory
        if firings_cnt+numel(fired)>size(firings,1)
            firings=[firings; zeros(10*numel(Tseq),2)];
        end
        firings(firings_cnt+1:firings_cnt+numel(fired),1:2)=[Tseq(t)+0*fired, fired];  
        firings_cnt=firings_cnt+numel(fired);
    end
    % store values for plotting, only last Sim.plot_length interval
    if Tseq(t)>Sim.t_End-Sim.plot_length
        vi=vi+1;
        % membrane potential
        PlotData.vplot(1:nplot,vi)=v(ind_plot); 
        % input currents
        PlotData.iEplot(1:nplot,vi)=J(ind_plot,1:N_e)*F.f(1:N_e);
        PlotData.iIplot(1:nplot,vi)=J(ind_plot,N_e+1:N_e+N_i)*F.f(N_e+1:N_e+N_i);
        PlotData.iExtplot(1:nplot,vi)=muRun(ind_plot,1);
    end
end

fprintf('--- End of trial...\n');
toc
%---------------------------------------------------------------------------
if ~any(any(firings))
    fprintf('\n --- NO SPIKES GENERATED... \n');
else
    % find last spike in firings
    IndexEnd=find(firings(:,2)==0,1)-1;
    if isempty(IndexEnd)
        IndexEnd=size(firings,1);
    end
    all_firings=firings(1:IndexEnd,[1 2]);
end


function F=synaptic_trace(Tau,dt,N_e,N_i)

    F=struct();
    tau_sE=Tau.tausyn_e; % exc synaptic time (fall time)
    tau_sI=Tau.tausyn_i; % inh synaptic time (fall time)
    fexp=[repmat(exp(-dt/tau_sE),N_e,1); repmat(exp(-dt/tau_sI),N_i,1)]; % Multiplicative step (fp)
    fSpike=[repmat((1/tau_sE),N_e,1); repmat((1/tau_sI),N_i,1)]; % add to fp with a spike
    f=zeros(N_e+N_i,1);
    F.fexp=fexp;
    F.fSpike=fSpike;
    F.f=f;

end

function F=syn_evolve(F,fired)


% update synaptic filter
    F.f=F.fexp.*F.f;
    if ~isempty(fired)      
    % update of synaptic filter with spikes
        F.f(fired)=F.f(fired)+F.fSpike(fired);
    end
end


end