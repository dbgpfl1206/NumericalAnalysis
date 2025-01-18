m=101;
v=1.0;
dx=1.0;
dt=0.02;
nmax=2500;
c(1)=1;

for i=1:m
    x(i)=(i-1)*dx;
end
%Exact Solution
for i=1:m
    d=v*dx/pe;
    yy1=(x(i)-v*dt*nmax)/(2*sqrt(d*dt*nmax));
    yy2=(x(i)+v*dt*nmax)/(2*sqrt(d*dt*nmax));
    ce(i)=c(1)/2*(erfc(yy1)+exp(v*x(i)/d)*erfc(yy2));
end
%print plots
subplot(1,3,1)
x=1:dx:m; 
y1=cc1(x); y2=cc2(x); y3=cc3(x); 
y4=cc4(x); y5=cc5(x); y6=cc6(x); 
y7=ce(x);
plot(x,y1,'-xb',x,y2,'-or',x,y3,'-+g',x,y4,'-dc',x,y5,'-vm',x,y6,'-py',x,y7,'-*k'); 
axis([1 100 0 1.2])
xlabel('Distance'),ylabel('Concentration')
title('pe=10')
legend('alpha=0, theta=0','alpha=0, theta=0.5','alpha=0, theta=1.0','alpha=1.0, theta=0','alpha=1.0, theta=0.5','alpha=1.0, theta=1.0','Exact solution')

%Calculation of direct RMSE and MSE
sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; sum6=0;
a1=0; a2=0; a3=0; a4=0; a5=0; a6=0;
for i=1:m
    a1=(ce(i)-cc1(i)).^2;
    sum1=sum1+a1;
end
for i=1:m
    a2=(ce(i)-cc2(i)).^2;
    sum2=sum2+a2;
end
for i=1:m
    a3=(ce(i)-cc3(i)).^2;
    sum3=sum3+a3;
end
for i=1:m
    a4=(ce(i)-cc4(i)).^2;
    sum4=sum4+a4;
end
for i=1:m
    a5=(ce(i)-cc5(i)).^2;
    sum5=sum5+a5;
end
for i=1:m
    a6=(ce(i)-cc6(i)).^2;
    sum6=sum6+a6;
end

%plot MSE
mse1=sum1/m; mse2=sum2/m; mse3=sum3/m; mse4=sum4/m; mse5=sum5/m; mse6=sum6/m;
mse=[mse1, mse2, mse3, mse4, mse5, mse6];
subplot(1,3,2) 
bar(mse)
title('pe=10'),ylabel('MSE')

%plot RMSE
rmse1=sqrt(sum1/m); rmse2=sqrt(sum2/m);
rmse3=sqrt(sum3/m); rmse4=sqrt(sum4/m);
rmse5=sqrt(sum5/m); rmse6=sqrt(sum6/m);
rmse=[rmse1, rmse2, rmse3,rmse4,rmse5,rmse6];
subplot(1,3,3)
bar(rmse)
title('pe=10'),ylabel('RMSE')