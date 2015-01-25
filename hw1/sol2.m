%% plot histograms

dist_name=cell(1,3);

dist_name{1}='gaussian';
dist_name{2}='1_-1';
dist_name{3}='cauchy';
maxIter=10000;
n=30;m=50;
for dist=1:3
    disp(['----' dist_name{dist}]);
    miuA=zeros(1,maxIter);
    ubA=zeros(1,maxIter);
    for iter=1:maxIter
        switch dist
            case 1
                A=randn(n,m);
            case 2
                A=(rand(n,m)>0.5)*2-1;
            case 3
                A=trnd(1,n,m);
        end
        A=normc(A);
        gramA=abs(A'*A);
        gramA(eye(m)==1)=0;
        miuA(iter)=max(gramA(:));
        ubA(iter)=1/2*(1+1/miuA(iter));
    end
    
    [h,x]=hist(miuA);
    h=h/maxIter;
    figure
    bar(x,h);
    title('histogram of  \mu(A)','FontSize',15,'FontWeight','Bold')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/' char(dist_name(dist)) '_h_mu.eps'] ,'epsc');

    [h,x]=hist(ubA);
    h=h/maxIter;
    figure
    bar(x,h);
    title('histogram of the upper bound','FontSize',15,'FontWeight','Bold')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/' char(dist_name(dist)) '_h_ub.eps'] ,'epsc');

end