%% SV: function description

clear all;
%function [est_theta] = SV(start, seg)
% [nfram, nbin, nchan]  = size(X); 
% X is a multichannel spectrum.



addpath ../../utils;

udir='testdata/';
edir='test_enh/';
regul=1e-3; % MVDR regularization factor
uname='10c_3m_90c';
part = strsplit(uname, 'c_');
incidence = str2num(part{1});

[x,fs]=audioread([udir uname '.wav']);

start=7;
seg=1;
sec = (start*fs+1:(start+seg)*fs);
x = x(sec,:);
nsampl=length(x);
nw=400;
inc = nw/2;
win = hanning(nw);
X = stft_multi(x.',win,inc);



% CGMM need [nfram, nchan, nbin]
X = permute(X, [2,3,1]);

[nfram, nchan,nbin]  = size(X); 

FS = [0:nbin-1]/nw*fs;
Xm = repmat(mean(X,1),nfram,1,1);
Xn = X - Xm;
% [L,~,Q]=CGMM_EM(Xn,2);
% plot(Q, '-*');
% save L L;
% 
% nv = permute(sum(L,1),[2,3,1]); %[nchan, nfram ]
% % 
% v = 2; % decision for the part of noise propobility.
% 
% steer=zeros(nchan,nbin);
% 
% for f = 1:nbin
% 	Rx = X(:,:,f)'*X(:,:,f)/nfram;
% 	Rn = zeros(nchan,nchan);
% 	for j = 1:nfram
% 		Rn = Rn + X(j,:,f)'*X(j,:,f)*L(j,v,f);
% 	end
% 	Rn  = Rn/nv(v,f);
% 	[V,D] = eig(Rx - Rn);
% 	[~,I] = max(diag(D));
% 	steer(:,f)=V(:,I);
% end 
% 
% save steer.mat steer;

% load steer.mat steer;
theta= incidence/180*pi;

array = [ -1,0;0,1;1,0;0,-1]*0.35;

interval = 0.35; % meter
c = 340; %meter/second.
TDOA = [-cos(theta), sin(theta),cos(theta),-sin(theta)]'*interval/c;
steer_ideal=sqrt(1/nchan)*exp(-2*1i*pi*TDOA*[0:nbin-1]/nw*fs); % steering vector

% D = 0:180;
% G = zeros(nbin,length(D));
% for i = 1:nbin
%     G(i,:)= beampattern(steer_ideal(:,i).',array,[i-1]*fs/nw,D);
% end

%mesh( D,[0:nbin-1]/nw*fs, G);
%plot(G(31,:))
phase = unwrap(angle(steer'));
figure,plot(FS, phase, '*-');
% 
% phase = angle(steer).';
% figure(2);plot(FS,phase(:,2), '*-r');
% phase = unwrap(phase,[],1);
% figure(3);plot(FS,phase(:,2),'*-r');
% %phase = angle(S{i}).';
% Wc = -2*pi*(0:nbin-1)*fs/nw;
% Ph = phase/P;
% in = Wc'\Ph;
% est_theta = atan2(in(1),in(2))/pi*180;
% disp(est_theta);
% 
% % 
% for i = 2:nbin
%     phase = angle(steer(:,i)).';
%     Wc = -2*pi*(i-1)*fs/nw;
%     Ph = phase/P;
%     in = Wc'\Ph;
%     sc(i) = atan2(in(1),in(2))/pi*180;
%     %disp(sc(i));
% end

disp('this is the end');

