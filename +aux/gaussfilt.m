% function [tout z] = gaussfilt(t,x,binwidth,causal_filter)
%
% filter 1D vector x with a Gaussian filter of std = binwidth (half width is
% half of that) by convolution of x with a gaussian. It uses 'conv'.
% Based on
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/272556
%
% INPUT:
% t = time vector
% x = data
% binwidth = twice the halfwidth (= std of gaussian filter). NOTE: must be
% integer > 1; width of gaussian kernel for filtering (same units as t);
% set it to 1 to prevent filtering  
% causal_filter = flag: if ==1, causal filter used (i.e., =0 at negative
% times).
%
% NOTE: a fictitious continuation of the t and x values to the left and
% right of the data vector are added to prevent bad filtering at the
% boundaries.   
%
% GLC, Nov 7, 2013

function  [tout z] = gaussfilt(t,x,binwidth,causal_filter)


if size(t,1)==1 t=t'; end;
if size(x,1)==1 x=x'; end;
exlen=binwidth;
x1=repmat(x(1),exlen,1);
x2=repmat(x(end),exlen,1);
x0=[x1 ; x ; x2];
tbin=diff(t(1:2));
t1=fliplr(t(1):-tbin:t(1)-(exlen-1)*tbin);
t2=fliplr(t(end):tbin:t(end)+(exlen-1)*tbin);
t0=[t1' ; t ; t2'];

% Construct blurring window:
windowWidth = int16(binwidth);
halfWidth = windowWidth / 2;
gaussFilter = gausswin(binwidth);
if causal_filter==1; gaussFilter(1:floor(max(size(gaussFilter))/2))=0; end;
if causal_filter==-1; gaussFilter(floor(max(size(gaussFilter))/2)+1:end)=0; end;
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

% alpha filter -- IN PROGRESS (13 nov 2013) -- do not trust
% bw=round(binwidth/2);
% alphac=1/bw/bw; y=alphac^2*(1:2*bw/alphac).*exp(-alphac*(1:2*bw/alphac));
% y=[ zeros(length(y),1)' y];
% masky=1:floor(length(y)/binwidth):length(y);
% y=y(masky);
% gaussFilter=y/sum(y);

% Do the blur:
z = conv(x0, gaussFilter);
z = z(halfWidth+exlen:end-halfWidth-exlen+1);
% tout = t0(halfWidth+exlen:end-halfWidth-exlen+1);
tout=exlen*tbin+t0(1:length(z));


% plot:
% figure(1); clf;
% hold on;
% plot(t,x,'k');
% plot(t0,z, 'r-', 'linewidth', 3);
% hold off;