[B,A] = ellip(4,1,20,0.5); % design lowpass filter
B = B .* (0.95).^[1:length(B)]; % contract zeros by 0.95
[H,w] = freqz(B,A);      % frequency response
theta = angle(H);        % phase response
thetauw = unwrap(theta); % unwrapped phase response