PlotParam=[];
for itrial=1:ntrials
    figure(3); clf;
    filename=fullfile(hmmdir,['Plot_trial' num2str(itrial) '.pdf']);
    PlotParam=struct('win',win_train(itrial,:),'colors',colors,'fntsz',8,'gnunits',gnunits);
    DATA=struct('win',win_train(itrial,:),'seq',hmm_postfit(itrial).sequence,'pstates',...
        hmm_results(itrial).pStates,'rates',hmm_results(itrial).rates);
    DATA.Spikes(1:gnunits)=spikes(itrial,1:gnunits);
%     figure(1); clf;
    hmm.fun_HMMRasterplot(DATA,HmmParam,PlotParam);
    saveas(gcf,filename,'pdf');
end