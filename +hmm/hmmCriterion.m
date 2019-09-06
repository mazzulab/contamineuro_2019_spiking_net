function outCrit=hmmCriterion(options,LL,varargin)

% input struct array dim [numel(VarStates),NumRuns], fields
% [tmp,epm,LLtrain]

if nargin>2
    temp_hmm_all_data=varargin{1};
end
SELECTION=options.SELECTION;
if strcmp(SELECTION,'min')
    [y_min,s_min]=hmm.funCriterion(LL,'min');
    stats=struct('LL_mean',LL);
elseif strcmp(SELECTION,'elbow')
    [y_min,s_min]=hmm.funCriterion(LL','elbow');
    stats=struct('LL_mean',LL);
elseif strcmp(SELECTION,'simMatch')
    matchValueTot=cell(1,size(temp_hmm_all_data,2));
    for i_st=1:size(temp_hmm_all_data,2)
        epm=arrayfun(@(x)x.epm,temp_hmm_all_data(:,i_st),'uniformoutput',0); % rows=different i.c.; cols=VarStates
        for i_run1=1:numel(epm)
            for i_run2=i_run1+1:numel(epm)
                % for all pairs of runs, match states by highest corr
                % is distribution of similarities bimodal? if yes, stop!
                rij=corr(epm{i_run1}',epm{i_run2}');
                matchSet=[1:size(rij,1); stableMatching(rij,rij)']';
                matchValue=rij(sub2ind(size(rij),matchSet(:,1)));
                matchValueTot{i_st}=[matchValueTot{i_st}; matchValue];
            end
        end
    end
    ave=cell2mat(cellfun(@(x)mean(x),matchValueTot,'uniformoutput',0));
    [y_min,s_min]=funCriterion(ave,'elbow');
    stats=struct('LL_mean',ave,'simAll',matchValueTot);
end
outCrit=struct('y_min',y_min,'s_min',s_min,'stats',stats);


















