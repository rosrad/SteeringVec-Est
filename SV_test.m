
clear all;
V = [4, 7, 11, 17.2, 20, 23]; %speech valid start point
D = [0, 5, 9, 12, 15.5 ]; %noise present start point

n=5; %
ValidSteer=zeros(n,2);
for j = 1:2
    for i=1:n
        ValidSteer(i,j) = SV(D(i),j);
    end
end
