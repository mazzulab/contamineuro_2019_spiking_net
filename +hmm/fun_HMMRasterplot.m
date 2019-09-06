% rasterplot with HMM
% OUTPUT h=figure handle

function fun_HMMRasterplot(DATA,HmmParam,PlotParam)

% VARIABLES
%     figure(1); 
% clf;

win=DATA.win; % xlimits on spikes plot
seq=DATA.seq;
pstates=DATA.pstates;
lambda=DATA.rates;
% PARAMETERS
BinSize=HmmParam.BinSize;
% colors=PlotParam.colors;
colshade=PlotParam.colors;
% colors=[colors; colors];
gnunits=PlotParam.gnunits;
fntsz=PlotParam.fntsz;
% legends
NumStates=size(pstates,1);
LegendUnits=cell(1,max(NumStates,gnunits));
for n=1:max(NumStates,gnunits)
    LegendUnits{n}=sprintf('%s',num2str(n));
end
AdjustX=0;


time_bins_pstates=DATA.win(1):BinSize:DATA.win(2)-BinSize; % time bins corresponding to pstates, already includes AdjustX offset
if numel(time_bins_pstates)<size(pstates,2)
    time_bins_pstates=[time_bins_pstates time_bins_pstates(end)+BinSize];
end
% remove offset from pstates
% hfig=figure(1); clf;
% hfig=figure(1); clf;
% set(gcf,'Visible','off');
%
PlotCols=1;
if ~isempty(seq)
    % admissible states
    State_spont=unique(seq(4,:),'Stable');
    PlotCols=numel(State_spont)+1;
end
subplot(2,PlotCols,[1 PlotCols])
%-------------
% STATES
%-------------
% seq
if ~isempty(seq)
    h=[];
    a=zeros(1,size(pstates,2));
    b=ones(1,size(pstates,2));
        for i=1:size(seq,2)
            ind=(time_bins_pstates>=seq(1,i) & time_bins_pstates<seq(2,i));
            [~,~]=aux.jbfill(time_bins_pstates(ind),...
                b(ind),a(ind),colshade(seq(4,i),1:3),0,0,0.1);
%                 b(ind),a(ind),colors(seq(4,i),:),0,0,0.1);
            hold on;
        end
%     end
    plot_cnt=0;
    for st=State_spont
        plot_cnt=plot_cnt+1;
%         size(time_bins_pstates)
%         size(pstates)
        h(plot_cnt)=plot(time_bins_pstates,pstates(st,:),'Color',colshade(st,:),'LineWidth',1.5);
        hold on
    end
    leg_letters=LegendUnits(State_spont);
    hold on
    xlab='Time [s]';
    ylab='Neurons';
    set(gca,'Yticklabel','none');
    ylim([0 1]);
    aux.figset(gca,xlab,ylab,'',fntsz)
end
% states that are not admissible
if isempty(seq)
    State_spont=[];
end
YLim=[0 1];

hold on;
% overlay spike train
Spikes=DATA.Spikes;

SpikeBar=0.3;
for unit=1:gnunits
    spikes_temp=[]; spks=[]; start_spk=[]; end_spk=[];
    start_spk=0+SpikeBar:1:gnunits-1+SpikeBar;
    end_spk=1-SpikeBar:1:gnunits-SpikeBar;
    start_spk=start_spk/gnunits;
    end_spk=end_spk/gnunits;
    spikes_temp=Spikes(unit).spk;

    spks=spikes_temp(spikes_temp>win(1) & spikes_temp<win(2));
    if ~isempty(spks) && numel(spks)~=2 
        line([spks spks]+AdjustX,[start_spk(unit) end_spk(unit)],'Color','k','LineWidth',1);
    elseif ~isempty(spks) && numel(spks)==2
        line([spks spks; spks spks]+AdjustX,[start_spk(unit) end_spk(unit)],'Color','k','LineWidth',1);
    end
    hold on
    set(gca,'Ytick',mean([start_spk; end_spk]),'YTickLabel',LegendUnits(1:gnunits));
end
xlim(win);
newPosition = [1-1/(numel(State_spont)+1)+0.02 0.25 0.05 0.05];
% legend(h,eventlabel(indlabel),'fontsize',8,'location',newPosition);

% PLOT FIRING RATES FOR admissible states only
if ~isempty(seq)
    lambda_plot=[]; x_shade=[];
    XLim=[0 1.1*max(max(lambda(State_spont(:),:)))];
    if any(any(lambda))
        for st=1:numel(State_spont)
%             if lfp_mode
%                 subplot(3,PlotCols,2*PlotCols+st)
%             elseif ~lfp_mode
                subplot(2,PlotCols,PlotCols+st)
%             end                
            [~,msg]=aux.jbfill(XLim,[gnunits+0.5 gnunits+0.5],zeros(1,2),colshade(State_spont(st),:),0.2,0,0.5);
            hold on
            y=barh(1:gnunits,lambda(State_spont(st),:));
            grey=[0.2 0.2 0.2];
            set(y,'facecolor',grey);
            % color
            if  0
                ylab='Neurons';
                xlab='Firing rate [spks/s]';
            else
                xlab='';
                ylab='';
            end
            aux.figset(gca,xlab,ylab,'',fntsz);
            hold on;
            text(XLim(2)*0.7,gnunits*0.9,leg_letters(st),...
                'LineWidth',1.5,'FontSize',10,'Color','k');
            if st==1
                set(gca,'Ytick',1:gnunits,'YTickLabel',LegendUnits(1:gnunits));
            else
                set(gca,'yticklabel',[],'xticklabel',[]);

            end
            xlim(XLim);
            ylim([0.5 gnunits+0.5]);
        end
    end
end
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 6 8 4]);
% p=mtit(plotTitle);
% filename=strrep(filename,'pdf','ai');
% set(gcf, 'renderer', 'painters');
% print( gcf, '-painters', filename, '-dill');


                

                