% Creates parameter files from choice of parameter set specified in Opt
% and saves it in .../DATA/Params.mat

function create_params_EI(paramsfile)
%----------
% NETWORKS
%----------
% Start and end of trial (in units of seconds)
Sim.t_Start=-1;
Sim.t_End=1;
Sim.dt_step=0.0001; % integration step (s)

%
%------------------
% NETWORK OPTIONS
%------------------
Network.clusters={'EE','EI','IE','II'}; % all weights are clustered
% Network.clust='hom'; % homogeneous EE cluster size
Network.clust='het'; % heterogeneous EE cluster size
Network.clust_std=0.01; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
Network.clustEI='hom'; % EI clusters: homogeneous='hom', heterogeneous='het'
Network.clustIE='hom'; % IE clusters: homogeneous='hom', heterogeneous='het'
Network.clustII='hom'; % II clusters: homogeneous='hom', heterogeneous='het'
% Network.clust_syn='';
N=2000; % network size
N_e = N*4/5; % exc neurons
N_i = N/5; % inh neurons
Scale=(1000/N)^(1/2);
% % global spontaneous firing rates (neeed to fix thresholds)
% ni_e = 2; % 3.6;   %6.6   % AB97: 3.0
% ni_i = 5; % 5.2;   %8.2   % AB97: 4.2
%------------------
% TIME CONSTANTS
%------------------
tau_arp = .005;  % refractory period
tau_i = .02;     % inh membrane time
tau_e = .02;	 % exc membrane time
tausyn_e=0.005;  % exc synaptic time
tausyn_i=0.005;  % inh synaptic time
%------------------
% SYNAPTIC WEIGHTS
%------------------
% syn weights are drawn from a gaussian distribution with std delta and
% mean Jab
Jee = Scale*0.02;
Jii = 1.*Scale*0.12; %     Jii = Scale*0.06;
Jie = 2*Scale*0.010; %     Jei 3 Scale*0.045;
Jei = 3.*Scale*0.02;
%------------------
% CLUSTER PARAMETERS
%------------------
% delta = 0.01;
delta=0.;% SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Jplus = 13; % EE intra-cluster potentiation factor
Network.factorEI = 10; % EI intra-cluster potentiation factor
Network.factorIE = 8; % IE intra-cluster potentiation factor
Network.factorII = 5; % II intra-cluster potentiation factor
%
bgr=0.1; % fraction of background excitatory neurons (unclustered)
Network.bgrI=0.1; % fraction of background neurons
Ncluster=100; % average # neurons per cluster
p = round(N_e*(1-bgr)/Ncluster); % # of clusters
f = (1-bgr)/p;
Network.fI = (1-Network.bgrI)/p;       % 0.09
% gam = 0.5;% % parameter related to inter-cluster depression (see function aux.SynWeights)
% rho = 2.75;
%------------------
% THRESHOLDS
%------------------
theta_e=1.42824; 
theta_i=0.74342;
% theta_e=1; % exc threshold potential
% theta_i=1; % inh threshold potentials
% reset potentials
He = 0;%
Hi = 0;%
%------------------
% CONNECTIVITY PARAMETERS
%------------------
Cee = N_e*0.2; % # presynaptic neurons
Cie = N_e*0.5; %
Cii = N_i*0.5; %
Cei = N_i*0.5; %
%------------------
% EXTERNAL BIAS
%------------------
% external input parameters, eg: external current given by mu_e_ext=Cext*Jee_ext*ni_ext
Cext = (N_e)*0.2; % # presynaptic external neurons
Jie_ext=0.8*Scale*0.0915;% external input synaptic strengths
Jee_ext=0.8*Scale*0.1027; %
% EXTERNAL CURRENT
% default external currents
ni_ext = 5; % 7;
mu_E0=Cext*Jee_ext*ni_ext;
mu_I0=Cext*Jie_ext*ni_ext;
% random ext current, different in each trial
Mu=[mu_E0*(ones(N_e,1)+(0.1/2)*(2*rand([N_e,1])-1)); ...
    mu_I0*(ones(N_i,1)+(0.05/2)*(2*rand([N_i,1])-1))];     % bias


%----------------
% DEFAULT STIMULI
%----------------
% STIMULUS
Stimulus.input='Const'; % constant external current
scnt=0;
% TASTE (specific stimulus)
scnt=scnt+1;
feat(scnt).name='US'; % unconditioned stimulus (taste)
feat(scnt).interval=[0 Sim.t_End]; % stimulus interval
gain=0.2; % stimulus value at 1 s
feat(scnt).gain=gain;
feat(scnt).profile=@(t)t; % time course of stimulus, eg a linear ramp
feat(scnt).selectivity='mixed';
feat(scnt).selective=rand(1,p)<0.5; % US selective clusters
feat(scnt).connectivity=0.5; % fraction of selective neurons within a selective cluster
% ANTICIPATORY CUE
scnt=scnt+1;
feat(scnt).name='CSgauss'; % conditioned stimulus (cue) with "quenched" noise
feat(scnt).interval=[-0.5 Sim.t_End]; % stimulus interval
gain=0.2; % SD of quenched noise across neurons
feat(scnt).gain=gain;
tau_cue=[0.5,1]; % rise and decay time of double exp cue time course
feat(scnt).profile=@(t)(1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1))); % double exp profile time course
feat(scnt).selectivity='exc'; % targets exc neurons only
feat(scnt).selective=ones(1,p); % CS targets all clusters
feat(scnt).connectivity=0.50; % fraction of neurons targeted within each cluster

Stimulus.feat=feat;


%------------------------------------------------------------------------
% PLOT PARAMETERS: ------------------------------------------------
%------------------------------------------------------------------------/
% PLOTS
Sim.Plotf=0;
Sim.plot_length=Sim.t_End-Sim.t_Start; % length of plot intervals
% indices of ensemble units to store
exc=randperm(N_e);
inh=N_e+1:N_e+N_i;
inh=inh(randperm(numel(inh)));
Sim.ind_p=[exc(1) inh(1)]; % choosing neuron index for membrane potential plot (one E and one I)
Sim.weights_save='off'; % save weight matrix: 'Yes'
extra='';
%
%

save(paramsfile,'ni_ext','tau_arp','tau_i','tau_e','theta_e',...
    'theta_i','delta','f','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext',...
    'Jplus','He','Hi','N_e',...
    'N_i','Cee','Cie','Cei','Cii','Cext','p','Sim','Network','Stimulus',...
    'tausyn_e','tausyn_i','extra','Mu','paramsfile');
% save(paramsfile,'ni_e','ni_i','ni_ext','tau_arp','tau_i','tau_e','theta_e',...
%     'theta_i','delta','f','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext',...
%     'Jplus','He','Hi','N_e',...
%     'N_i','Cee','Cie','Cei','Cii','Cext','p','Sim','Network','Stimulus',...
%     'tausyn_e','tausyn_i','extra','Mu','paramsfile');

fprintf('Network parameters saved in %s\n',paramsfile);