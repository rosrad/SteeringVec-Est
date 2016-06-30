close all; clear;
addpath ../../voice;
addpath ../../utils;
chime_data='G:\corpus\CHiME\chime3';
upath=[chime_data, '/data/audio/16kHz/isolated/']; % path to segmented utterances
epath=[chime_data, '/data/audio/16kHz/CGMM/']; % path to enhanced utterances
cpath=[chime_data, '/data/audio/16kHz/embedded/']; % path to continuous recordings
bpath=[chime_data, '/data/audio/16kHz/backgrounds/']; % path to noise backgrounds
apath=[chime_data, '/data/annotations/']; % path to JSON annotations
%nchan=6;

% Define hyper-parameters
pow_thresh=-20; % threshold in dB below which a microphone is considered to fail
wlen = 1024; % STFT window length
regul=1e-3; % MVDR regularization factor
cmin=6400; % minimum context duration (400 ms)
cmax=12800; % maximum context duration (800 ms)

sets={'dt05'};
modes={'real'};
chanlist=[1 3:6];
nchan=length(chanlist);
for set_ind=1:length(sets),
    set=sets{set_ind};
    for mode_ind=1:length(modes),
        mode=modes{mode_ind};
        
        % Read annotations
        mat=json2mat([apath set '_' mode '.json']);
        real_mat=json2mat([apath set '_real.json']);
            
        for utt_ind=1:length(mat),
                           
            if ~strcmpi(mat{utt_ind}.environment, 'bus')
                continue;
            end
            
            
            udir=[upath set '_' lower(mat{utt_ind}.environment) '_' mode '/'];
            edir=[epath set '_' lower(mat{utt_ind}.environment) '_' mode '/'];
            if ~exist(edir,'dir'),
                system(['mkdir -p ' edir]);
            end
            uname=[mat{utt_ind}.speaker '_' mat{utt_ind}.wsj_name '_' mat{utt_ind}.environment];
            
            % Load WAV files
            xsize=audioinfo([udir uname '.CH' int2str(chanlist(1)) '.wav']);
            nsampl=xsize.TotalSamples;
            x=zeros(nsampl,nchan);
            for c=1:nchan,
                [x(:,c),fs]=audioread([udir uname '.CH' int2str(chanlist(c)) '.wav']);
            end
            
            % Check microphone failure
            xpow=sum(x.^2,1);
            xpow=10*log10(xpow/max(xpow));
            fail=(xpow<=pow_thresh);
            
            % Load context (up to 5 s immediately preceding the utterance)
            if strcmp(mode,'real'),
                cname=mat{utt_ind}.wavfile;
                cbeg=round(mat{utt_ind}.start*16000)-cmax;
                cend=round(mat{utt_ind}.start*16000)-1;
                for utt_ind_over=1:length(mat),
                    cend_over=round(mat{utt_ind_over}.end*16000);
                    if strcmp(mat{utt_ind_over}.wavfile,cname) && (cend_over >= cbeg) && (cend_over < cend),
                        cbeg=cend_over+1;
                    end
                end
                cbeg=min(cbeg,cend-cmin);
                n=zeros(cend-cbeg+1,nchan);
                for c=1:nchan,
                    n(:,c)=audioread([cpath cname '.CH' int2str(chanlist(c)) '.wav'],[cbeg cend]);
                end
            elseif strcmp(set,'tr05'),
                cname=mat{utt_ind}.noise_wavfile;
                cbeg=round(mat{utt_ind}.noise_start*16000)-cmax;
                cend=round(mat{utt_ind}.noise_start*16000)-1;
                n=zeros(cend-cbeg+1,nchan);
                for c=1:nchan,
                    n(:,c)=audioread([bpath cname '.CH' int2str(chanlist(c)) '.wav'],[cbeg cend]);
                end
            else
                cname=mat{utt_ind}.noise_wavfile;
                cbeg=round(mat{utt_ind}.noise_start*16000)-cmax;
                cend=round(mat{utt_ind}.noise_start*16000)-1;
                for utt_ind_over=1:length(real_mat),
                    cend_over=round(real_mat{utt_ind_over}.end*16000);
                    if strcmp(mat{utt_ind_over}.wavfile,cname) && (cend_over >= cbeg) && (cend_over < cend),
                        cbeg=cend_over+1;
                    end
                end
                cbeg=min(cbeg,cend-cmin);
                n=zeros(cend-cbeg+1,nchan);
                for c=1:nchan,
                    n(:,c)=audioread([cpath cname '.CH' int2str(chanlist(c)) '.wav'],[cbeg cend]);
                end
            end
            
            % STFT
            nw=400;
            %win = ones(nw,1);
            win = hanning(nw);
            X = stft_multi(x.',win,nw/4);
            X = permute(X,[2,3,1]);
            %X=stft_multi(x.',wlen);
            %X = permute(X,[2,3,1]);
            [nfram,nchan,nbin] = size(X);    
            % Compute noise covariance matrix
            N = stft_multi(n.',win,nw/4);
            N = permute(N,[2,3,1]);
            Ncov=zeros(nchan,nchan,nbin);
            for f=1:nbin,
                for i=1:size(N,1),
                    Ntf=permute(N(i,:,f),[2 1 3]);
                    Ncov(:,:,f)=Ncov(:,:,f)+Ntf*Ntf';
                end
                Ncov(:,:,f)=Ncov(:,:,f)/size(N,1);
            end

%             Xcov = zeros(nchan,nchan,nbin);
%             for f=1:nbin
%                 Xcov(:,:,f) = X(:,:,f)'*X(:,:,f)/size(X,1); 
%             end
            %norm_X = X - repmat(mean(X,1),size(X,1),1);
            [L,R,Q]=CGMM_EM(X,2);
            nv = permute(sum(L,1),[2,3,1]);
            v = 2;
            steer=zeros(nchan,nbin);
            for f = 1:nbin
                Rx = X(:,:,f)'*X(:,:,f)/nfram;
                Rn = zeros(nchan,nchan);
                for j = 1:nfram
                    Rn = Rn + X(j,:,f)'*X(j,:,f)*L(j,v,f);
                end
                Rn  = Rn/nv(v,f);
                [V,D] = eig(Rx - Rn);
                [~,I] = max(diag(D));
                steer(:,f)=V(:,I);
            end 
            %load( 'v2.new.mat', 'steer');
            % Localize and track the speaker
            %[~,TDOA]=localize(permute(X,[3,1,2]),chanlist);
            
            % MVDR beamforming
            Xspec=permute(mean(abs(X).^2,1),[2 3 1]);
            Y=zeros(nbin,nfram);
            
            %Xcov = Ncov(~fail,~fail,f)+regul*diag(Xspec(~fail,f))
            for f=1:nbin,
                Df =steer(:,f);
                for t=1:nfram,
                    Xtf=permute(X(t,:,f),[2 3 1]);
                    %Df=sqrt(1/nchan)*exp(-2*1i*pi*(f-1)/wlen*fs*TDOA(:,t)); % steering vector
                    Y(f,t)=Df(~fail)'/(Ncov(~fail,~fail,f)+regul*diag(Xspec(~fail,f)))*Xtf(~fail)/(Df(~fail)'/(Ncov(~fail,~fail,f)+regul*diag(Xspec(~fail,f)))*Df(~fail));
                    %Y(f,t)=Df(~fail)'/Xcov(~fail,~fail,f)*Xtf(~fail)/(Df(~fail)'/Xcov(~fail,~fail,f)*Df(~fail));
                end
            end
            %y = istft_multi(Y', win,nw/4,nsampl);
            y=istft_multi(Y,nsampl,win,nw/4).';
            % Write WAV file
            y=y/max(abs(y));
            audiowrite([edir uname 'v2.Ncov' '.wav'],y,fs);
            return;
        end
    end
end
 
