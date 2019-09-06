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


function [J, params]=fun_SynWeights_EI(paramsfile)

% LOAD PARAMETERS

params=load(paramsfile);
aux.v2struct(params);
Network=params.Network;
aux.v2struct(Network);
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
gam=1/(2-f*(Q+1));%
Jminus = 1.-gam*f*(Jplus-1.);
params.Jminus = Jminus;

% INHIBITORY CLUSTERS
if any(strcmp(clusters,'EI'))
    %     JminusEI = 1.-gam*fI*(JplusEI-1.);
    JplusEI = 1/(1/p+(1-1/p)/factorEI);
    JminusEI = JplusEI/factorEI;
    params.JminusEI = JminusEI;
    jei_out=-JminusEI*Jei; % intra-cluster
    jei_in=-JplusEI*Jei; % inter-cluster
    fprintf('JplusEI=%0.03g, JminusEI=%0.03g\n',JplusEI,JminusEI);
end
if any(strcmp(clusters,'IE'))
    %     JminusIE = 1.-gam*fI*(JplusIE-1.);
    JplusIE = 1/(1/p+(1-1/p)/factorIE);
    JminusIE = JplusIE/factorIE;
    params.JminusIE = JminusIE;
    jie_out=JminusIE*Jie; % intra-cluster
    jie_in=JplusIE*Jie; % inter-cluster
    fprintf('JplusIE=%0.03g, JminusIE=%0.03g\n',JplusIE,JminusIE);
end
if any(strcmp(clusters,'II'))
    JplusII=factorII;
    JminusII = 1.-gam*fI*(JplusII-1.);
    params.JminusII = JminusII;
    jii_out=-JminusII*Jii; % intra-cluster
    jii_in=-JplusII*Jii; % inter-cluster
    fprintf('JplusII=%0.03g, JminusII=%0.03g\n',JplusII,JminusII);
end

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
peeout=pee;
peein=pee;
peiout=pei;
peiin=pei;
pieout=pie;
piein=pie;
piiout=pii;
piiin=pii;
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
JEI=(jei*(ones(N_e,N_i)+deltaEI*randn(N_e,N_i))).*(rand([N_e,N_i])<pei);
JIE=(jie*(ones(N_i,N_e)+deltaIE*randn(N_i,N_e))).*(rand([N_i,N_e])<pie);
JII=(jii*(ones(N_i)+delta*randn(N_i,N_i))).*(rand([N_i,N_i])<pii);
clustermatrix=eye(Q);
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


%-------------
% EI CLUSTERS
%-------------
% check #clusters and coding level are consistent
if any([strcmp(clusters,'EI'), strcmp(clusters,'EI'), strcmp(clusters,'II')])
    NcUnits=round(fI*N_i);    %  number of Exc units per cluster
    fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
    Numbg=round(N_i*(1-fI*p)); % number of background (i.e. non-selective) Exc units
    fprintf('  --- fraction of bg Inh units: %0.03g',Numbg/N_i);
    switch Network.clustEI
        case 'hom'
            popsizeI=repmat(NcUnits,Q,1);
        case 'het'
            Nc=[];
            clust_std=Network.clust_std;
            while (sum(Nc)-(N_i-Numbg))~=0 || any(Nc<0)
                Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
            end
            popsizeI=Nc; % array of cluster sizes
            if any(sum(popsizeI)-(N_i-Numbg))
                fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
            end
    end
    cusumNcI=[0 cumsum(popsizeI)'];
    
    %-------------
    % EI weights
    %-------------
    if any(strcmp(Network.clusters,'EI'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustEI,'het') || strcmp(Network.clustEI,'hom')
            % clustered units: inter-cluster weights
            JEI(1:cusumNcE(Q+1),1:cusumNcI(Q+1))=...
                (jei_out*(ones(cusumNcE(Q+1),cusumNcI(Q+1))+deltaEI*randn(cusumNcE(Q+1),...
                cusumNcI(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcI(Q+1)])<peiout); % inter-cluster weights
            for clu=2:Q+1 % intra-cluster higher weights
                JEI(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcI(clu-1):cusumNcI(clu))=...
                    (jei_in*(ones(popsize(clu-1),popsizeI(clu-1))+deltaEI*randn(popsize(clu-1),popsizeI(clu-1)))).*...
                    (rand([popsize(clu-1),popsizeI(clu-1)])<peiin);
            end
        end
    end
    
    
    %-------------
    % EI CLUSTERS
    %-------------
    %-------------
    % IE weights
    %-------------
    if any(strcmp(Network.clusters,'IE'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustIE,'het') || strcmp(Network.clustIE,'hom')
            % clustered units: inter-cluster weights
            JIE(1:cusumNcI(Q+1),1:cusumNcE(Q+1))=...
                (jie_out*(ones(cusumNcI(Q+1),cusumNcE(Q+1))+deltaIE*randn(cusumNcI(Q+1),...
                cusumNcE(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcE(Q+1)])<pieout); % inter-cluster weights
            for clu=2:Q+1 % intra-cluster higher weights
                JIE(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcE(clu-1):cusumNcE(clu))=...
                    (jie_in*(ones(popsizeI(clu-1),popsize(clu-1))+deltaIE*randn(popsizeI(clu-1),popsize(clu-1)))).*...
                    (rand([popsizeI(clu-1),popsize(clu-1)])<piein);
            end
        end
    end
    
    %-------------
    % II CLUSTERS
    %-------------
    
    %-------------
    % II weights
    %-------------
    if any(strcmp(Network.clusters,'II'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustII,'het') || strcmp(Network.clustII,'hom')
            % clustered units: inter-cluster weights
            JII(1:cusumNcI(Q+1),1:cusumNcI(Q+1))=...
                (jii_out*(ones(cusumNcI(Q+1),cusumNcI(Q+1))+deltaII*randn(cusumNcI(Q+1),...
                cusumNcI(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcI(Q+1)])<piiout); % inter-cluster weights
            for clu=2:Q+1 % intra-cluster higher weights
                JII(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcI(clu-1):cusumNcI(clu))=...
                    (jii_in*(ones(popsizeI(clu-1),popsizeI(clu-1))+deltaII*randn(popsizeI(clu-1),popsizeI(clu-1)))).*...
                    (rand([popsizeI(clu-1),popsizeI(clu-1)])<piiin);
            end
        end
    end
    params.popsizeI=popsizeI;
end






JEI(JEI>0)=0;
JIE(JIE<0)=0;
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
