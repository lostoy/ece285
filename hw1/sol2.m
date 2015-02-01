%% plot histograms

dist_name=cell(1,3);

dist_name{1}='gaussian';
dist_name{2}='1_-1';
dist_name{3}='cauchy';
maxIter=1000;
n=30;m=50;
miuA=zeros(1,maxIter);
ubA=zeros(1,maxIter);
for iter=1:maxIter
    A=randn(n,m);
    A=normc(A);
    gramA=abs(A'*A);
    gramA(eye(m)==1)=0;
    miuA(iter)=max(gramA(:));
    ubA(iter)=1/2*(1+1/miuA(iter));
end

[h,x]=hist(miuA,100);
h=h/maxIter/(x(2)-x(1));

figure
% appoximation result from 

%Cai, Tony, Jianqing Fan, and Tiefeng Jiang. 
%"Distributions of angles in random packing on spheres." 
%The Journal of Machine Learning Research 14.1 (2013): 1837-1864.

a=m^(2/(n-1));
K=1/4/sqrt(pi)*gamma(n/2)/gamma((n+1)/2);
tx=0:0.01:1-0.01;
pdf_x=K*exp(-K*(a*acos(tx)).^(n-1)).*(n-1).*(a*acos(tx)).^(n-2).*a./sqrt(1-(tx).^2)+...
    K*exp(-K*(a*acos(-tx)).^(n-1)).*(n-1).*(a*acos(-tx)).^(n-2).*a./sqrt(1-(-tx).^2);
plot(tx,pdf_x,'r--',x,h,'b.','LineWidth',2)

title(['histogram of  \mu(A),n=30,m=' num2str(m)],'FontSize',15,'FontWeight','Bold')
set(gca,'FontSize',15,'FontWeight','Bold');
legend('asymptotic pdf','histogram')
saveas(gca, ['./eps/' 'm_' num2str(m) '_h_mu.eps'] ,'epsc');

[h,x]=hist(ubA,100);
h=h/maxIter/(x(2)-x(1));
figure
plot(x,h,'LineWidth',2);
title('histogram of the upper bound','FontSize',15,'FontWeight','Bold')
set(gca,'FontSize',15,'FontWeight','Bold');
saveas(gca, ['./eps/' 'm_' num2str(m) '_h_ub.eps'] ,'epsc');



%%
