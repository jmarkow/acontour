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


% visualization parameters

font_size=12;
disp_band=[0 12e3];
clipping=-10;
imscale=15;
timescale_list = 0.5:0.2:2.2; % time scales in milliseconds for this analysis
sono_choice=8;

% Load the sample sound
[signal, fs] = wavread('finchdoublet.wav');

tic
[consensus_power,f,t,auditory_contours,sonograms]=acontour(signal,fs,'timescale_list',timescale_list);
toc

% to pass parameters use parameter/value pairs e.g.
%
% [consensus_power,f,t,auditory_contours,sonograms]=acontour(signal,fs,'len',30);
%
% would change the window length to 30 ms


%% Visualization

figure
h(1) = subplot(2,1,1);
imagesc(t,f,max(log(consensus_power+imscale),clipping));
ylim([disp_band])
ylabel('Frequency [Hz]','FontSize', font_size);
axis xy;
xlabel('Time [s]','FontSize', font_size);
title('Consensus contour representation','FontSize', font_size);
colormap(hot)

h(2) = subplot(2,1,2);
imagesc(t,f,max(log(abs(sonograms{sono_choice})+imscale),clipping));
ylim([disp_band]);
titlev=sprintf('Sonogram computed with timescale %6.3g ms',timescale_list(sono_choice));
title(titlev,'FontSize', font_size);
ylabel('Frequency [Hz]','FontSize', font_size);
set(gca,'YDir','normal');
xlabel('Time [s]','FontSize', font_size);
colormap(hot)

linkaxes(h,'x')




