function temp_sequence=fun_HMM_postfit_sequence(pstates_rate,HmmParam,win_coding)

BinSize=HmmParam.BinSize;
MinDur=round(HmmParam.MinDur/(BinSize));
MinP=HmmParam.MinP;
numstates=size(pstates_rate,1);
temp_sequence=[]; 
for st_cnt=1:numstates
    % 0's where this holds, 1's otherwise
    x=~(pstates_rate(st_cnt,:)>MinP);
    % find all start and end bins and duration of each interval
    dsig = diff([1 x 1]);
    startIndex = find(dsig < 0);
    endIndex = find(dsig > 0)-1;
    duration = endIndex-startIndex+1;
    % keep only events with duration>MinDur
    stringIndex = find(duration >= MinDur);
    % revision: keep 1st states spilling over from
    % pre-delivery, but then ignore their initial
    % transition in HMM3_rates later
    % create sequence data
    % only states lasting more than MinDur ms get included in the sequence 
    tot_time=[]; % collect time interval when admissible states are pstates>HmmParam.MinP;
    if ~isempty(stringIndex)
        feat_temp=zeros(numel(stringIndex),5);
        for ind_st=stringIndex
            time=[]; temp_startIndex=[]; temp_endIndex=[];
            temp_startIndex=startIndex(ind_st);
            temp_endIndex=endIndex(ind_st);
            time=(temp_endIndex(1)-temp_startIndex(1)+1)*BinSize;
            % collect all start and end times for all states in temp_sequence
            % 1st row=start (s); 
            % 2nd row=end (s); 
            % 3rd row=duration
            % 4th row=state as in pstate
            temp_sequence=[temp_sequence [(temp_startIndex-1)*BinSize+win_coding(1);...
                temp_endIndex*BinSize+win_coding(1); time; st_cnt]];

        end
    end
end
if ~isempty(temp_sequence)
    [~, I]=sort(temp_sequence(1,:));
    temp_sequence=temp_sequence(:,I);
end
