
% This program demonstrates the recently developed time-frequency method to 
% represent the dynamically changing sound signal like bird song with 
% auditory contours selected at their own natural parameter sets, time-scale and angle.
% To reduce the computational load, we provide an approximated algorithm 
% for defining auditory object in time-frequency plane. This program 
% extracts auditory contours following well-developed image processing 
% methods for connecting neighboring pixels in 2D space.
% Although this program is coded to reduce the computational load, 
% we recommend to use the short sound signal(3-4 secs) unless if you have a large memory space.

% To run this program, signal processing, image processing and parallel computing toolbox are required.
% If your machine contains multiple cores of CPU, open a pool of MATLAB workers.
% To do that, type "matlabpool open" in the command window

% For more detail method, refer
% Y Lim, BG Shinn-Cunningham, and TJ Gardner, Stable Time-Frequency Contours for Sparse Signal Representation, EUSIPCO-2013

% Authors: Yoonseob Lim, Jeff Markowitz, and Timothy J. Gardner
% Copyright (C) Yoonseob Lim, Jeff Markowitz, and Timothy J. Gardner
% All rights reserved


tic

% Load the sample sound
[signal, fs] = wavread('finchdoublet.wav');

	
%% Parameters

FREQUENCYLIMIT = 12e3; %Upper Range of plot in Hz

%these parameter settings are more practical, for speed.
Nfft = 1024;
Nshift = 20;

angles = (pi/8:pi/8:pi) + pi/8; % compute contours in all these angles - recommended not to change.
TScale = 0.5:0.2:2.2; % time scales in milliseconds for this analysis
sonogramindex = 8; % which of the time-scales above should be used for the sonogram?
ARThreshold = 95; % keep only contours longer than this length percentile 98 or 99 recommended.

NtScale = length(TScale);
Nangle = length(angles);

tWin = -Nfft/2+0.5:Nfft/2-0.5;

%% Visualization

PIXELRANGE=floor((Nfft/2)*1000*FREQUENCYLIMIT/(fs/2));
FONTSIZE=15;

% compute the Time and Frequency vector for plotting

figure
h(1) = subplot(2,1,1);
imagesc(T,F,log(consensuspower+10))
ylim([0 FREQUENCYLIMIT])
ylabel('frequency [kHz]','FontSize', FONTSIZE);
set(gca,'YDir','normal');
xlabel('time [ms]','FontSize', FONTSIZE);
title('Consensus contour representation','FontSize', FONTSIZE);
colormap(hot)

h(2) = subplot(2,1,2);
imagesc(T,F,log(abs(sonogramrecord{sonogramindex})+10))
ylim([0 FREQUENCYLIMIT])
titlev=sprintf('Sonogram computed with timescale %6.3g ms',TScale(sonogramindex));
title(titlev,'FontSize', FONTSIZE);
ylabel('frequency [kHz]','FontSize', FONTSIZE);
set(gca,'YDir','normal');
xlabel('time [ms]','FontSize', FONTSIZE);
colormap(hot)

linkaxes(h,'x')
toc




