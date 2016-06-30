%D= 0:180;G= beampattern(ones(1,5)/5, ULA*0.2, 1000, D);figure,plot(D,G)

function G = beampattern(wop, array,fs, D)
m =size(array,1);
%fs =1000;
am = zeros(m,length(D));
%wop=ones(1,m)/m;
L = -2:2;
for k = 1:length(D)
    h = [cosd(D(k)), sind(D(k))].';
    coefficient = 2 * pi * fs/ 340;
    %tau1 = L*0.2*sind(90-D(k));
    tau = array * h;
    am(:,k) = exp(j * coefficient * tau);
end

G = 10*log10(abs(wop*am));
%A=abs(wop*am);
%G=20*log10(abs(wop*am));
%plot(A);
% mesh(degree,FS,A');
% 
% function beampattern()
% FS = 1:2000;
% D = 0:180;
% array = [ -1,0;0,1;1,0;0,-1];
% %array = [ -2,0;-1,0;0,0; -1,0; -2,0];
% 
% G = zeros(length(D),length(FS));
% parfor i =1:length(FS)
%     G(:,i) = fs_beam(array,FS(i), D);
% end
% mesh(FS,D,G);