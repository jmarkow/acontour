function [CONSENSUS,F,T,AUDITORY_CONTOUR,SONOGRAM]=acontour_core(SIGNAL,FS,varargin)
% This program demonstrates the recently developed time-frequency method to 
% represent the dynamically changing sound SIGNAL like bird song with 
% auditory contours selected at their own natural parameter sets, time-scale and angle.
% To reduce the computational load, we provide an approximated algorithm 
% for defining auditory object in time-frequency plane. This program 
% extracts auditory contours following well-developed image processing 
% methods for connecting neighboring pixels in 2D space.
% Although this program is coded to reduce the computational load, 
% we recommend to use the short sound SIGNAL(3-4 secs) unless if you have a large memory space.

% To run this program, SIGNAL processing, image processing and parallel computing toolbox are required.
% If your machine contains multiple cores of CPU, open a pool of MATLAB workers.
% To do that, type "matlabpool open" in the command window

% For more detail method, refer
% Y Lim, BG Shinn-Cunningham, and TJ Gardner, Stable Time-Frequency Contours for Sparse Signal Representation, EUSIPCO-2013

% Authors: Yoonseob Lim, Jeff Markowitz, and Timothy J. Gardner
% Copyright (C) Yoonseob Lim, Jeff Markowitz, and Timothy J. Gardner
% All rights reserved

if nargin<1 | isempty(SIGNAL)
	error('Need signal to continue');
end

if nargin<2 |  isempty(FS)
	error('Need sampling rate to continue');
end



%% PARAMETERS (PASSED AS PARAMETER/VALUE PAIRS) 

%these parameter settings are more practical, for speed.

len = 23.2; % window length (in ms)
overlap = 22.8; % window overlap (in ms)

angle_list = (pi/8:pi/8:pi) + pi/8; % compute contours in all these angle_list - recommended not to change.
timescale_list = 0.5:0.2:2.2; % time scales in milliseconds for this analysis
clength_threshold = 95; % keep only contours longer than this length percentile 98 or 99 recommended.
norm_amp=1;
filtering=[];
pow_weight=1;
zeropad=0;

%% END USER PARAMETERS

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'len'
			len=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'angle_list'
			angle_list=varargin{i+1};
		case 'timescale_list'
			timescale_list=varargin{i+1};
		case 'clength_threshold'
			clength_threshold=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'norm_amp'
			norm_amp=varargin{i+1};
		case 'pow_weight'
			pow_weight=varargin{i+1};
		case 'zeropad'
			zeropad=varargin{i+1};
	end
end

len=round((len/1e3)*FS);
overlap=round((overlap/1e3)*FS);

if zeropad==0
	zeropad=round(len/2);
end

if ~isempty(zeropad)
	SIGNAL=[zeros(zeropad,1);SIGNAL(:);zeros(zeropad,1)];
	disp(['Zero pad: ' num2str(zeropad/FS) ' S']);
end

if norm_amp
	disp('Normalizing signal amplitude');
	SIGNAL=SIGNAL./max(abs(SIGNAL));
end

if ~isempty(filtering)
	disp(['Filtering the signal (' num2str(filtering) ' Hz) high pass']);
	[b,a]=ellip(5,.2,40,[filtering]/(FS/2),'high');
	SIGNAL=filtfilt(b,a,SIGNAL);
end


% convert window and overlap into samples

ntimescale_lists = length(timescale_list);
nangle_list = length(angle_list);

t_win = -len/2+0.5:len/2-0.5;

%% Auditory contour
% Here we first extract all the contours at each parameter and select only 
% longer contours. The reason why we select longer contours is because those
% contours capture most part of the SIGNAL structure. We also found that
% longer contours are structurally stable even with a small perturbations
% on time-scale or angle_list. This means that this part of the code can be useful for
% the analysis of other kinds of time-series such as EEG or LFP if
% the analysis is to find the consistent structures in the SIGNAL over many
% trials.

disp('Computing contours...');

%AUDITORY_CONTOUR = cell(1,ntimescale_lists);
for sigmacount = 1:ntimescale_lists %use the matlab parallel computing toolbox if multiple processors are available.

	timescale = timescale_list(sigmacount);
	timescale = (timescale/1000)*FS;
	w = exp(-(t_win/timescale).^2); % gauss window
	dw = w.*(t_win/(timescale^2))*-2; % deriv gauss window

	q=spectrogram(SIGNAL,w,overlap,[],FS)+eps;
	q2=spectrogram(SIGNAL,dw,overlap,[],FS)+eps;
	dx = (q2./q)/(2*pi); %displacement according to the remapping algorithm
	SONOGRAM{sigmacount} = q;

	for angle_variable = 1:nangle_list

		theta = angle_list(angle_variable);

		% a quick approximation to the contour detection
		s = -1*(imag(dx*exp(1j*theta))<0)+(imag(dx*exp(1j*theta))>0);
		[gx, gy] = gradient(s);
		BW = ((-gx*cos(theta+pi/2)+gy*sin(theta+pi/2))>.001);        

		cc = bwconncomp(BW); %build a strucure containing a separate entry for each contour (connected component of BW).
		cc_pix = regionprops(cc,'Area'); %this is the length of each contour
		weightv = [];        
		for i = 1:length(cc_pix),
			weightv(i)=cc_pix(i).Area;
		end

		% Select only longer contours
		long_contours = find(weightv>=prctile(weightv,clength_threshold));

		% Construct 2d image that contains contour indices
		contour_image = zeros(size(q));        
		power_image = zeros(size(q));
		
		for ik = 1:length(long_contours)
			
			ind = cc.PixelIdxList(long_contours(ik));
			
			% integrate the pwr	

			pwr=sum(abs(q(ind{1})));
			contour_image(ind{1}) = 1;
			power_image(ind{1}) = pwr;

		end
		
		AUDITORY_CONTOUR(sigmacount,angle_variable).bin_img = sparse(contour_image);
		AUDITORY_CONTOUR(sigmacount,angle_variable).pwr_img = sparse(power_image);
		
	end

end

disp('Computing consensus...');

%% Consensus process
% Superimpose all contours from all angle_list and time scales to produce a single image
% Through this process, contours with stable structures will be highlited.

CONSENSUS = zeros(size(SONOGRAM{1}));

for sigmacount = 1:ntimescale_lists-1

	consensus = zeros(size(CONSENSUS));   
	for angle_variable=1:nangle_list

		if angle_variable==1
			cv = AUDITORY_CONTOUR(sigmacount,1).bin_img + ...
				AUDITORY_CONTOUR(sigmacount+1,angle_variable).bin_img + ...
			       	AUDITORY_CONTOUR(sigmacount,nangle_list).bin_img;
			consensus = consensus + (cv>1);
		else         
			cv = AUDITORY_CONTOUR(sigmacount,angle_variable).bin_img + ...
				AUDITORY_CONTOUR(sigmacount+1,angle_variable).bin_img + ...
			       	AUDITORY_CONTOUR(sigmacount,angle_variable-1).bin_img;
			consensus = consensus + (cv>1); %keep only contour points that show some agreement across neighboring angle_list or time-scales.
		end
	end

	if pow_weight
		CONSENSUS = CONSENSUS + consensus.*abs(SONOGRAM{sigmacount}); % multiply consenus score by sonogram power for final image
	else
		CONSENSUS = CONSENSUS + consensus;
	end

end


% compute the Time and Frequency vector for plotting

F = linspace(0,FS/2,len/2+1);
sig_len=length(SIGNAL);
ncol = fix((sig_len-overlap)/(len-overlap));
colindex = 1 + (0:(ncol-1))*(len-overlap);
T = ((colindex-1)+((len)/2)')/FS; 

if ~isempty(zeropad)
	T=T-zeropad/FS;
end
