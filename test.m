close all; clear;
addpath ../../voice;

chime_data='G:\corpus\CHiME\chime3';
upath=[chime_data, '/data/audio/16kHz/isolated/']; % path to segmented utterances

file='G:\corpus\CHiME\chime3\data\audio\16kHz\enhanced\dt05_bus_real\F01_050C010R_BUS.wav';
[x,fs] =audioread(file);
% 
% 
%  OV=4;                               % overlap factor of 2 (4 is also often used)
%  INC=100;                             % set frame increment in samples
%  NW=INC*OV;                          % DFT window length
%  W=sqrt(hamming(NW,'periodic'));     % omit sqrt if OV=4
%  W=W/sqrt(sum(W(1:INC:NW).^2));      % normalize window
%  F=rfft(enframe(S,W,INC),NW,2);      % do STFT: one row per time frame, +ve frequencies only
%  
%  X=overlapadd(irfft(F,NW,2),W,INC);
 
% INC=20;       						% set frame increment in samples
% NW=INC*2;
%     % oversample by a factor of 2 (4 is also often used)
% S=cos((0:NW*7)*6*pi/NW);				% example input signal 	% sqrt hamming window of period NW
% [F,t,w]=enframe(S,NW,INC);               	% split into frames
% X=overlapadd(F,w,INC);  
nw=40;
x=cos((0:nw*7)*6*pi/nw);	

%W = hamming(nw,'periodic');
INC = nw/4;
n = size(x,1);
[frm,t,W] = enframe(x,nw,INC, 'z');
%F=rfft(frm,nw,2);      % do STFT: one row per time frame, +ve frequencies only
X=overlapadd(frm,W,INC);
X = X(1:n);

y = X-x;
s = y> 0.00000000001;
%X = stft_multi(x, hamming(nw),nw/4,'z');

%xx = istft_multi(X, hamming(nw),nw/4,size(x,1));