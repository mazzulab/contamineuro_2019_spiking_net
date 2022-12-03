% function fun_HMM(METHOD,spikes)
%
% Run HMM on one session
%
% Luca Mazzucato April 2016

function OUT=fun_HMM_modelSel(IN)

aux.v2struct(IN);
% with fields: 'METHOD','spikes','win_train','HmmParam'
[ntrials, gnunits]=size(Spikes);

% HMM PARAMETERS
temp_SkipSpikesSess=[];
HmmParam=struct();
HmmParam.AdjustT=0.; % interval to skip at trial start to avoid canonical choise of 1st state in matlab
HmmParam.BinSize=0.002;%0.005; % time step of Markov chain
HmmParam.MinDur=0.05;   % min duration of an admissible state (s) in HMM DECODING
HmmParam.MinP=0.8;      % pstate>MinP for an admissible state in HMM ADMISSIBLE STATES
HmmParam.NumSteps=10;%    %10 number of fits at fixed parameters to avoid non-convexity
HmmParam.NumRuns=50;%     % 50% % number of times we iterate hmmtrain over all trials
SELECTION='elbow';
switch METHOD
    case 'XVAL'
        NP=@(NumStates,gnunits,logT)0;
        HmmParam.NumSteps=1;
        HiddenStep=1; % increment in hidden states% RUNNING EACH SESSION SEPARATELY
        HiddenMin=4;
        HiddenMax=50; % max value allowed
    case 'BIC'
        NP=@(NumStates,gnunits,logT)(NumStates.*(NumStates-1)+NumStates.*gnunits+NumStates-1)*logT;
        fprintf('\n BIC...\n');
        HiddenStep=1; % increment in hidden states% RUNNING EACH SESSION SEPARATELY
        HiddenMin=2;
        HiddenMax=50; % max value allowed
    case 'AIC'
        NP=@(NumStates,gnunits,logT)(NumStates.*(NumStates-1)+NumStates.*gnunits+NumStates-1)*2;
        fprintf('\n AIC...\n');
        HiddenStep=1; % increment in hidden states% RUNNING EACH SESSION SEPARATELY
        HiddenMin=2;
        HiddenMax=50; % max value allowed
    case 'simMatch'
        NP=@(NumStates,gnunits,logT)0;
        fprintf('\n Similarity matching...\n');
        HiddenStep=1; % increment in hidden states% RUNNING EACH SESSION SEPARATELY
        HiddenMin=2;
        HiddenMax=50; % max value allowed
end
HmmParam.NP=NP;
optionsCrit=struct('SELECTION',SELECTION);
if any(strcmp(fieldnames(IN),'SELECTION'))
    optionsCrit=struct('SELECTION',IN.SELECTION);
end

% ADAPTIVE INCREASE OF HIDDEN STATES
% if the mean LL in the last two VarStates is decreasing, keep
% adding a new batch of hidden states
% if it increased, stop
%
LLtot=struct('m2LL',[]);
stop_LL=0;
HiddenTotal=[]; % store all numbers of hidden states
temp_hmm_all_data=[]; % only save in BIC-AIC
temp_all_bestfit=[]; % only save in BIC-AIC

s_old=1;
count_scan=0; % need at least 5 steps with same criterion stopping point to select it
outCrit=struct();
% if XVAL, generate training and test sets on top
myCluster = parcluster('local');
NumWorkers=myCluster.NumWorkers;
if strcmp(METHOD,'XVAL')
    K=20; % # of xval runs
    fL=0.5; % fraction of holdouts, take at most half of trial
    L=round(fL*ntrials); % don't let L be larger than # of trials
    K=min([K,ntrials]); % don't let K be larger than # of trials
    fprintf('leave-%d-out XVAL parameters:\n',L);
    fprintf('%d xval runs\n',K);
    %
    % indices for XVAL
    indxval=repmat(struct('train',[],'test',[]),1,K);
    for k=1:K
        a=randperm(ntrials);
        indxval(k).train=a(1:ntrials-L);
        indxval(k).test=a(ntrials-L+1:ntrials);
    end
    % adjust parfor batches to NumWorkers
    GoodSteps=max(2,floor((NumWorkers-3)/K));
else
    % adjust parfor batches to NumWorkers
    GoodSteps=max(2,floor((NumWorkers-3)/HmmParam.NumSteps));
end
HiddenBatch=HiddenMin:HiddenStep:HiddenMin+HiddenStep*(GoodSteps-1);
while stop_LL==0
    tic
    HiddenTotal=[HiddenTotal HiddenBatch];
    HmmParam.VarStates=HiddenBatch;
    switch METHOD
        case 'XVAL'
            [LLtotxval,hmm_all_data,temp_SkipSpikesSess]=hmm.fun_HMM_XVAL(Spikes,indxval,win_train,HmmParam);
            temp_hmm_all_data=[temp_hmm_all_data hmm_all_data];
            LLtot.m2LL=[LLtot.m2LL; -2*LLtotxval'];
        otherwise % BIC or AIC
            [LLtottemp,hmm_all_data,hmm_all_bestfit,temp_SkipSpikesSess]=hmm.fun_HMM_BIC_AIC(Spikes,win_train,HmmParam);
            temp_hmm_all_data=[temp_hmm_all_data hmm_all_data];
            temp_all_bestfit=[temp_all_bestfit hmm_all_bestfit];
            LLtot.m2LL=[LLtot.m2LL; LLtottemp];
    end
    toc
    % CRITERIONS:
    % 1) min; 2) elbow; 3) drop in similarity of best matching features across NumRuns
    if length(HiddenTotal)>4
        %         % start checking criterion only after slope drops
        %         if diff(Ydiff(i_sm-2:i_sm-1))<std(Ydiff(1:i_sm-1))/10 && startCrit==0; startCrit=1; end
        %             [~,s_shift] = funCriterion(Y(2:i_sm-1),opt.ModelSel);  % skip first entry in Y
        if strcmp(optionsCrit.SELECTION,'simMatch')
            outCrit = hmm.hmmCriterion(optionsCrit,LLtot.m2LL,temp_hmm_all_data);
        else
            outCrit = hmm.hmmCriterion(optionsCrit,LLtot.m2LL);
        end
        s_min=outCrit.s_min;
        if s_min==s_old
            count_scan=count_scan+numel(HiddenBatch);
            fprintf('found min -> count %d\n',count_scan);
        else
            count_scan=1;
        end
        s_old=s_min;
    end
    Last=HmmParam.VarStates(end);
    if count_scan>6
        fprintf('Best fit found, end of training...');
        break;
    else
        if Last<HiddenMax
            HiddenBatch=Last+HiddenStep:HiddenStep:Last+HiddenStep*GoodSteps;
            fprintf('LL still decreasing (-2LL=%0.03g), keep training up to %d states...\n',LLtot.m2LL(end),HiddenBatch(end));
        elseif Last>=HiddenMax
            stop_LL=1;
            fprintf('Reached upper limit %d on # of states, end of training...\n',Last);
        end
    end
end
StatesSelected=HiddenTotal(s_min);
fprintf('%d states...\n',StatesSelected);

fieldNames={'HiddenTotal',...
    'HmmParam','temp_hmm_all_data','temp_SkipSpikesSess','LLtot','StatesSelected','outCrit','fieldNames'};
OUT=aux.v2struct(fieldNames);
