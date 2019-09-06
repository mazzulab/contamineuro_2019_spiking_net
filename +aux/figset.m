% figset(hh,xlab,ylab,tt,fntsz,varargin)
%
% utility function which provides my favorite figure settings
% 
% - hh: current figure handle (e.g., gca)
% - xlab: (string) xlabel
% - ylab: (string) ylabel
% - tt:   (string) title
% - fntsz: (numeric) font size (typically between 15 and 22)
% - varargin{1}: (string) zlabel in 3D plots
% 
% USAGE: figset(gca,xlab,ylab,tt,fntsz,varargin);
%
% GLC, 20 Nov 2007 -- added optional color parameter Aug 2012

function figset(hh,xlab,ylab,tt,fntsz,varargin)
col='none'; 
xlabel(xlab,'fontsize',fntsz);
ylabel(ylab,'fontsize',fntsz);
title(tt,'fontsize',fntsz);
set(hh,'fontsize',fntsz,'color',col);
box off;
if ~isempty(varargin) % 2D plot
    zlab=varargin{1};
    zlabel(zlab,'fontsize',fntsz);
end
