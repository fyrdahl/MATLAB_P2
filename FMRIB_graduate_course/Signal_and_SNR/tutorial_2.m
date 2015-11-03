%	Part 2 - Signal Simulator
%	==============================================================
%%	Initial Simulator Plot
T1			=	1200;
TR			=	1000;
flipAngle	=	90;
tMax		=	5000;
simSignal;

%%	Flip angle 65
flipAngle 	=	65;
simSignal;

%%	Question 2.3
TR	=	20;
T1	=	1200;
flipAngle	=	10.5;
simSignal;

%%	Ernst Angle
simSignalvFlip;

%	Part 3 - SNR
%	==============================================================
%%	Interactive ROI SNR calculation
load structural;
SNR_calc_with_ROIs(structural);

%%	Histogram Definition
load masks;
tissue_roi		=	find(tissue_mask);
background_roi	=	find(background_mask);
tissue_vals		=	structural(tissue_roi);
background_vals	=	structural(background_roi);

S 	=	mean(tissue_vals);
N	=	std(tissue_vals);

tissue_bins		=	linspace(S-3*N, S+3*N, 50);
background_bins	=	linspace(-3*N, 3*N, 50);

figure();
subplot(1,2,1);
hist(tissue_vals, tissue_bins);
xlim([S-2*N, S+2*N]);

subplot(1,2,2);
hist(background_vals, background_bins);
xlim([-2*N, 2*N]);

%%	Interactive ROI SNR calculation
SNR_calc_with_ROIs(structural);
SNR_calc_with_ROIs(structural, 1000);

%%	EPI Average Data
load epi_data;
SNR_calc_with_ROIs(epi_ave1);
SNR_calc_with_ROIs(epi_ave9);
