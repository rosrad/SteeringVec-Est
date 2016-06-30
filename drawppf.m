function y=drawpp(m,wop,f,myarray)
%figure,plot(e)
thetas=[0:180];
tm=thetas*pi/180;

am= zeros(m,181);
for iDegree = 0:180
     h = [cosd(iDegree),sind(iDegree)].';
     coefficient = 2 * pi * f/ 340;
     tau = myarray * h;
     am(:,iDegree+1) = exp(j * coefficient * tau);
end


% am=exp(j*2 * pi * f * tau);
A=abs(wop'*am);  %阵列响应
A=A/max(A);
figure,polar(tm,A)
A=10*log10(A);  %对数图
%hold on,title('归一化阵列响应幅值极坐标图，八阵元，信噪比20db')
figure,plot(thetas,A);
%hold on,title('八阵元，信噪比20db')
%hold on,xlabel('入射角/度')
%hold on,ylabel('归一化 A=10*log10(A);')
grid on 
axis([0 180 -35 0]);
%hold on,plot(-45,-35:0.1:0,'r');
%hold on,plot(30,-35:0.1:0,'r');
%hold on,plot(0,-35:0.1:0,'r');
%hold on,plot(60,-35:0.1:0,'r');
hold on,plot(20,-35:0.1:0,'r');
hold on,plot(60,-35:0.1:0,'r');
