function res=funHMM(DATAIN)

RUN=1;
Plotf=1;
METHOD=DATAIN.METHOD;
if strcmp(METHOD,'XVAL')
    SELECTION='elbow'; % 'min';%
else
    SELECTION='min';
end
% if any(strcmp(fieldnames(DATAIN),'SELECTION'))
%     SELECTION=DATAIN.SELECTION;
% end
% if any(strcmp(fieldnames(DATAIN),'RUN'))
%     RUN=DATAIN.RUN;
% end
spikes=DATAIN.spikes;
gnunits=size(spikes,2);
win_train = DATAIN.win;
% run HMM
IN=struct('METHOD',METHOD,'Spikes',spikes,'win_train',win_train,'SELECTION',SELECTION);
OUT=hmm.fun_HMM_modelSel(IN);
%
aux.v2struct(OUT);
% PICK BEST FIT
BestStateInd=find(StatesSelected==HiddenTotal);
HmmParam.VarStates=StatesSelected;
colors=aux.distinguishable_colors(max(HiddenTotal(BestStateInd),4));
% POSTFIT ANALYSIS
% if XVAL, rerun training and decoding on all data
% if ~XVAL, pick best fit
hmm_results=[];
hmm_data=temp_hmm_all_data;
if strcmp(METHOD,'XVAL')
    HmmParam.NumSteps=10;
    [sequence, SkipSpikesSess]=hmm.fun_HMM_binning(spikes,HmmParam,win_train);
    temp_hmm_all_data0=hmm.fun_HMM_training(sequence,gnunits,HmmParam);
    % -> This is needed to avoid getting stuck in local minima (HMM is non-convex)
    tempLL=cell2mat(arrayfun(@(x)x.LLtrain,temp_hmm_all_data0(:,1),'uniformoutput',false));
    [~,ind_step]=min(tempLL); % find index of initial cond. with highest LL
    hmm_bestfit=temp_hmm_all_data0(ind_step); % epm fit
else
    SkipSpikesSess=OUT.temp_SkipSpikesSess;
    tempLL=cell2mat(arrayfun(@(x)x.LLtrain,temp_hmm_all_data(:,BestStateInd),'uniformoutput',false));
    [~,ind_step]=min(tempLL); % find index of initial cond. with highest LL
    hmm_bestfit=temp_hmm_all_data(ind_step,BestStateInd);
end
% CREATE SEQUENCES
hmm_results=hmm.fun_HMM_decoding(spikes,hmm_bestfit,HmmParam,win_train);
% HMM ADMISSIBLE STATES
hmm_postfit=hmm.fun_HMM_postfit(spikes,hmm_results,HmmParam,win_train);

% reinitialize in each session
hmm_multispikes=SkipSpikesSess;
HmmParam.VarStates=HiddenTotal;
res=aux.v2struct({'fieldNames','hmm_data','hmm_bestfit','hmm_results','hmm_postfit','hmm_multispikes','HmmParam','win_train','colors','BestStateInd','LLtot'});
if any(strcmp(fieldnames(DATAIN),'filesave'))
    filename=[DATAIN.filesave '_' METHOD '.mat'];
    save(filename,'res');
    fprintf('\n --- hmm data saved in %s. \n',filename);
end
%--------
% PLOTS
%--------
%%
figure(1); clf; h=[]; LLtemp=[];
for i=1:numel(HiddenTotal); LLtemp=[LLtemp,[hmm_data(:,i).LLtrain]']; end
SD=nanstd(LLtemp)/sqrt(size(LLtemp,1));
[~,msg]=aux.jbfill(HiddenTotal,LLtot.m2LL'+SD,LLtot.m2LL'-SD,[0.8,0.8,0.8],0,0.1); hold on;
h=plot(HiddenTotal,LLtot.m2LL(:,1)','color','k','linewidth',2); hold on;
plot(HiddenTotal(BestStateInd),LLtot.m2LL(BestStateInd,1),'o','color','r','markersize',10); hold on;
aux.figset(gca,'hidden states',sprintf('%s',METHOD),'',15);
filename=[DATAIN.filesave '_' METHOD '.pdf'];
saveas(gcf,filename,'pdf');
if 0
    % plot how many skipped spikes in binned sequences
    figure(1); clf;
    h=SpikesSkippedFun(SkipSpikesSess,gnunits);
    filename=[files.hmm_plot '_' an '_' sess 'SkippedSpikes' num2str(HmmParam.BinSize) 'msBin.pdf'];
    saveas(gcf,filename,'pdf');
end
