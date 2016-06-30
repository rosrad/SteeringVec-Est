close all; clear;
d = 2;
k = 3;
n = 500;
[X,label] = mixGaussRnd(d,k,n);

figure(1);
plotClass(X,label);

m = floor(n/2);
X1 = X(:,1:m);
X2 = X(:,(m+1):end);
% train
[z1,model,llh] = mixGaussEm(X1,k);
figure(2);
plot(llh);
figure(3);
plotClass(X1,z1);
% predict
z2 = mixGaussPred(X2,model);
figure(4);
plotClass(X2,z2);