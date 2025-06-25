
%% This function aims to convert single side bispectrum to double side

% Author: JL

% Different to ordinary methods one should not directly mulitply all
% positive components by 2 because the spectrum is not entirely symmetrical
% around 0

%% Functions part

function [faxis,bspec_single] = bispect_sing(waxis,bspec,Fs)

% Input: 
% 'waxis' and 'bspec' from bi-spectrum analysis
% Fs: sampling rate

% Output:
% the frequency axis and output single-sided bi-spectrum

idx_poszero = (waxis>=0);
idx_pos = (waxis>0);
idx_neg = (waxis<0 & waxis>-0.5);
idx_zero = (waxis==0);

faxis = waxis(idx_poszero).*Fs;

bspec_single = NaN(length(idx_poszero),length(idx_poszero));

% Extract 3 other domains
bspec_xneg_ypos = abs(bspec(idx_neg,idx_pos));
bspec_xpos_yneg = abs(bspec(idx_pos,idx_neg));
bspec_xneg_yneg = abs(bspec(idx_neg,idx_neg));

bspec_xzero_yneg = abs(bspec(idx_zero,idx_neg));
bspec_xneg_yzero = abs(bspec(idx_neg,idx_zero));

% Add up to main domain
bspec_poszero = abs(bspec(idx_poszero,idx_poszero));
% Add DC, flip relevant axes
bspec_poszero(2:end,1) = bspec_poszero(2:end,1) + flip(bspec_xneg_yzero);
bspec_poszero(1,2:end) = bspec_poszero(1,2:end) + flip(bspec_xzero_yneg);

% Add other frequency
%bspec_poszero(2:end,2:end) = bspec_poszero(2:end,2:end) + flip(bspec_xpos_yneg,1);
%bspec_poszero(2:end,2:end) = bspec_poszero(2:end,2:end) + flip(bspec_xneg_ypos,2);
bspec_poszero(2:end,2:end) = bspec_poszero(2:end,2:end) + flip(flip(bspec_xneg_yneg,1),2);



bspec_single = bspec_poszero;


