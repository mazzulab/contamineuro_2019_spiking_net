% function PlotHistStairs
%
% plot histogram of a distribution using stairs
% INPUT: DATA=raw data
%         NBins=number of histogram bins
%         colors='k','r',...
%         varargin{1}=OPTION: if OPTION='log', plot log scale, otherwise
%                       set OPTION=[];
%         varargin{2}=XLim: first and last bins
% OUTPUT: h=plot handle
%
% LM November 2015

function h=PlotHistStairs(DATA,NBins,OPTIONS)
% %
% colors=colors(1);
% NBins=15;
% varargin{1}='log';
% varargin{2}=XLim;
colors='k'; % default color
if any(strcmp(fieldnames(OPTIONS),'colors'))
    colors=OPTIONS.colors;
end
XLim=[];
if any(strcmp(fieldnames(OPTIONS),'XLim'))
    XLim=OPTIONS.XLim;
end
linestyle='-';
if any(strcmp(fieldnames(OPTIONS),'linestyle'))
    linestyle=OPTIONS.linestyle;
end

% histogram corr coeff
h=[];
% plot distribution
MIN=0.011;
DATA=DATA(~isnan(DATA));
Xmax=max(DATA);
Xmin=min(DATA);
if ~isempty(XLim)
    x=linspace(XLim(1),XLim(2)*1.01,NBins+1);
else
    x=linspace(Xmin,Xmax*1.01,NBins+1);
end
f=histc(DATA,x);
f=f'/(diff(x(1:2))*numel(DATA));
f=[0 f];
x=[x(1) x];
h=stairs(x,f,'linewidth',3);
hold on;
set(h,'color',colors,'linestyle',linestyle);
