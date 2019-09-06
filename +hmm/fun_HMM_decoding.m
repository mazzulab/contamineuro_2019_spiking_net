%---------------------
% HMM_decoding
%---------------------
% optional 3rd output: row vector with total LL of decoding for each
% VarStates

function hmm_all_results=fun_HMM_decoding(spikes,hmm_bestfit,HmmParam,win_train)
% hmm_bestfit=hmm_all_bestfit;
% skip this interval at the beginning of the spike train: needed to avoid
% the problem with hmmtrain in matlab that always starts with state 1 in
% the first bin.

% option for bernoulli or poisson (only affecting in firing rate estimation)
DISTR='bernoulli';
if any(strcmp(fieldnames(HmmParam),'distr'))
    DISTR=HmmParam.distr;
end

AdjustT=0.;
if any(strcmp(fieldnames(HmmParam),'AdjustT'))
    AdjustT=HmmParam.AdjustT; 
end
win_decode=win_train;
win_decode(:,1)=win_decode(:,1)-AdjustT;
% initialize variables
VarStates=HmmParam.VarStates;
BinSize=HmmParam.BinSize;
[ntrials,gnunits]=size(spikes);
hmm_all_results=repmat(struct('rates',[],'pStates',[],'Logpseq',[]),ntrials,numel(VarStates));
% START DECODING
% posterior state distributions
% compute the probability that at each step we are in state s, using the
% transition and emission matrices [esttr,estemis] estimated with hmmtrain above.
tr_cnt=0;
LL=10^10*ones(numel(VarStates),sum(ntrials)); % store log-likelihoods for comparison
seq=hmm.fun_HMM_binning(spikes,HmmParam,win_decode);

for trial=1:ntrials
    tr_cnt=tr_cnt+1;
    % pstates has: rows=states; columns=time bins
    % logpseq=log probability of the sequence given esttr and
    % estemis
    for st_cnt=1:numel(VarStates)
        NumStates=VarStates(st_cnt);
        esttr=hmm_bestfit(st_cnt).tpm;
        estemis=hmm_bestfit(st_cnt).epm;
        % hmmdecode return NaN if any entry in esttr or estemis is zero
        % fix it with this hack (incorporated in hmmdecode.m)
        MinValue=1e-100;
        esttr(esttr<MinValue)=MinValue;
        estemis(estemis<MinValue)=MinValue;
%         if strcmp(DISTR,'poisson')
%             [~,~,logpseq,pstates]=ph_em(seq(trial).data',esttr,estemis,1);
%             pstates=pstates';
%         else
            [pstates,logpseq]=hmmdecode(seq(trial).data,esttr,estemis);
%         end
        % FIT FIRING RATES IN EACH STATE
        % POST DELIVERY RATE FIT
        % remove first 0.2 s in pstates
        pstates_rate=pstates(:,round(AdjustT/BinSize)+1:end);
        v1=ones(size(pstates_rate,2),1);
        lambda=zeros(NumStates,gnunits);
        win=win_train(trial,:);
        for st=1:NumStates
            for unit=1:gnunits
                yt=zeros(1,size(pstates_rate,2));
                % turn spikes into yt
                if strcmp(DISTR,'poisson')
                    spks=hmm.Spikes2Seq(spikes(trial,unit),win,BinSize,DISTR); % get spike count
                    A=(pstates_rate(st,:)*spks')/(pstates_rate(st,:)*v1);
                    lambda(st,unit)=(A/BinSize);
                else
                    temp_spikes=spikes(trial,unit).spk;
                    spks=temp_spikes(temp_spikes>win(1) & temp_spikes<win(2))-win(1);
                    spks=ceil(spks/BinSize);
                    yt(spks)=1;
                    A=(pstates_rate(st,:)*yt')/(pstates_rate(st,:)*v1);
                    lambda(st,unit)=-(1/BinSize)*log(1-A);
                end
            end
        end
        hmm_all_results(trial,st_cnt).rates=lambda;
        hmm_all_results(trial,st_cnt).pStates=pstates_rate;
        hmm_all_results(trial,st_cnt).Logpseq=logpseq;
    end
end