% function [J, params]=SynWeights(params)
%
% OUTPUT
%       J                =synaptic matrix
%       OUT.popsize      =array of dim # pops, how many exc neurons in each
%                         pop
%       OUT.clustermatrix=matrix of dimension # pops x # clusters, each row shows
%                   which clusters belong to that pop (1s and 0s)
%
% Luca Mazzucato February 2016


function [J, params]=fun_SynWeights(paramsfile)

% LOAD PARAMETERS

params=load(paramsfile);
aux.v2struct(params);
Network=params.Network;
Sim=params.Sim;
%-----------------------
% PARAMETERS VALUES
%-----------------------
numfig=1;
Next=N_e; % external units
% CLUSTERS
Q=p; % number of clusters
%-----------------------
% SYNAPTIC WEIGHTS
%-----------------------
% WEIGHTS
% depression
gam = 0.5;% % parameter related to inter-cluster depression Jminus
Jminus = 1.-gam*f*(Jplus-1.);
params.Jminus = Jminus;
%
jee=Jee;
jee_out=Jminus*Jee; % intra-cluster potentiation
jee_in=Jplus*Jee; % inter-cluster depression
jei=-Jei;
jie=Jie;
jii=-Jii;
% connection probability
pee=Cee/N_e;
pie=Cie/N_e;
pei=Cei/N_i;
pii=Cii/N_i;
pext=Cext/Next;
%
fprintf('  --- Jplus=%0.03g, Jminus=%0.03g\n',Jplus,Jminus);

%----------------------------
% SYNAPTIC MATRIX
%----------------------------
% if delta>0, generate a distribution of synaptic weights with mean J and
% variance delta^2 J^2
peeout=pee;
peein=pee;
% check #clusters and coding level are consistent
NcUnits=round(f*N_e);    %  number of Exc units per cluster
Numbg=round(N_e*(1-f*p)); % number of background (i.e. non-selective) Exc units
switch Network.clust
    case 'hom'
        popsize=repmat(NcUnits,Q,1);
    case 'het'
        Nc=[];
        clust_std=Network.clust_std;
        while (sum(Nc)-(N_e-Numbg))~=0 || any(Nc<0)
            Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
        end
        popsize=Nc; % array of cluster sizes
        if any(sum(popsize)-(N_e-Numbg))
            fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
        end
end
cusumNcE=[0 cumsum(popsize)'];
% background units (if present), if not, override in next line
JEE=(jee*(ones(N_e)+delta*randn(N_e,N_e))).*(rand([N_e,N_e])<peeout);
if strcmp(Network.clust,'het') || strcmp(Network.clust,'hom')
    % clustered units: inter-cluster weights
    JEE(1:cusumNcE(Q+1),1:cusumNcE(Q+1))=...
        (jee_out*(ones(cusumNcE(Q+1))+delta*randn(cusumNcE(Q+1),...
        cusumNcE(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcE(Q+1)])<peeout); % inter-cluster weights
    for clu=2:Q+1 % intra-cluster higher weights
        JEE(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcE(clu-1):cusumNcE(clu))=...
            (jee_in*(ones(popsize(clu-1))+delta*randn(popsize(clu-1),popsize(clu-1)))).*...
            (rand([popsize(clu-1),popsize(clu-1)])<peein);
    end
end
clustermatrix=eye(Q);

JEI=(jei*(ones(N_e,N_i)+delta*randn(N_e,N_i))).*(rand([N_e,N_i])<pei);
JEI(JEI>0)=0;
JIE=(jie*(ones(N_i,N_e)+delta*randn(N_i,N_e))).*(rand([N_i,N_e])<pie);
JIE(JIE<0)=0;
JII=(jii*(ones(N_i)+delta*randn(N_i,N_i))).*(rand([N_i,N_i])<pii);
JII(JII>0)=0;
JEE(JEE<0)=0;
J=[JEE JEI; JIE JII];
J=J-diag(diag(J)); % eliminate self-couplings
fprintf('  --- New synaptic weights set: done...\n');
fprintf('      Overall: Jee=%g -- Jie=%g -- Jei=%g -- Jii=%g \n',jee,jie,jei,jii);
fprintf('      Var[J]=(Jx%0.03g)^2\n',delta);

params.popsize=popsize;
params.clustermatrix=clustermatrix;

if strcmp(Network.clust,'het')
    for i=1:Q
        a(i)=sum(popsize(clustermatrix(:,i)>0));
    end
    fprintf('  --- clusters size: mean=%0.03g neurons/cluster, sd/mean=%0.03g\n',mean(a),std(a)/mean(a));
end
fprintf('  --- fraction of bg Exc units: %0.03g\n',N_e/Numbg);

%-------------------
% PLOT weight matrix
%-------------------
figure(1); clf; 
subplot(2,1,1);
colormap(parula); %xy=J; fun_colormapLim;
imagesc(J);
aux.figset(gca,'neurons','neurons','weights',10);
colorbar;
subplot(2,1,2);
lambda=eig(J);
plot(real(lambda),imag(lambda),'.');
aux.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
saveas(gcf,fullfile('data','weights.pdf'),'pdf');

% spike thresholds for each population
theta=[params.theta_e, params.theta_i];
Theta=zeros(1,numel(params.popsize)+2);
Theta(1:numel(params.popsize)+1)=theta(1);
Theta(end) = theta(2);
params.Theta=Theta;
fprintf('Spike thresholds Theta calculated for each population and stored in params\n');
