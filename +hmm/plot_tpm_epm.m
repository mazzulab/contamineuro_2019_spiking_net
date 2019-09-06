tpm0=hmm_bestfit.tpm;
epm0=hmm_bestfit.epm/HmmParam.BinSize; epm0(:,end)=[]; %epm0=epm0';
col=colors;
% sort states from the longest lived
[B,I]=sort(diag(tpm0),'descend');
tpm=tpm0(I,I); epm=epm0(I,:); col=col(I,:);
% TPM
minExp=4;
tpm(tpm<10^(-minExp))=10^(-minExp);
figure(1); clf;
colormap('parula');
imagesc((log10(tpm)));
hC=colorbar;% set(hC,'ColorScale','log')
l = fliplr(-linspace(0,minExp,minExp+1)); % Tick mark positions
Lab=[]; Lab{1}=1;
for i_l=2:numel(l); Lab{i_l}=sprintf('10^{-%d}',i_l-1); end
ylabel(hC,'prob(i$\to$j)','FontSize',20,'interpreter','latex');
set(hC,'Ytick',l,'YTicklabel',fliplr(Lab));
aux.figset(gca,'State','State','tpm',20);
filename=fullfile(hmmdir,'Tpm.pdf');
saveas(gcf,filename,'pdf');
% EPM
figure(2); clf;
subplot(12,1,1);
for i_c=1:size(col,1)
    [~,~]=aux.jbfill([i_c-0.5,i_c+0.5],...
        ones(1,2),zeros(1,2),col(i_c,:),0,0,1); hold on
end
xlim([0.5,i_c+0.5]); title('raster colors'); set(gca,'ytick',[]);
subplot(12,1,[3 12]);
colormap('hot');
hi=imagesc(fliplr(epm)'); hC=colorbar;% set(hC,'ColorScale','log')
set(gca,'ytick',[1,size(epm,2)],'yticklabel',[size(epm,2),1]);
ylabel(hC,'Firing rate [spks/s]','FontSize',15);
aux.figset(gca,'State','Neuron','epm',15);
filename=fullfile(hmmdir,'Epm.pdf');
saveas(gcf,filename,'pdf');
