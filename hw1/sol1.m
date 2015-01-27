clc;
close all;
%% generating A and normalize
n=30;m=50;

%% eval random x for uniform, gaussian , {+1,-1}, cauchy
maxIter=1000;
maxSupp=10;
dist_name=cell(1,3);

dist_name{1}='uniform';
dist_name{2}='gaussian';
dist_name{3}='1_-1';
dist_name{4}='cauchy';

options.min_error=1e-4;
options.t=0.5;
runTimes=zeros(1,5);
for dist=1:4
    disp(['----' dist_name{dist}]);
    
    
    err1_MP=zeros(maxSupp,maxIter);
    err2_MP=zeros(maxSupp,maxIter);
    
    err1_OMP=zeros(maxSupp,maxIter);
    err2_OMP=zeros(maxSupp,maxIter);
    
    err1_LSOMP=zeros(maxSupp,maxIter);
    err2_LSOMP=zeros(maxSupp,maxIter);
    
    err1_WMP=zeros(maxSupp,maxIter);
    err2_WMP=zeros(maxSupp,maxIter);
    
    err1_ThMP=zeros(maxSupp,maxIter);
    err2_ThMP=zeros(maxSupp,maxIter);
    
    for supp_num=1:maxSupp
        disp(['-------- supp= ' num2str(supp_num) ' / ' num2str(maxSupp)])
        for iter=1:maxIter
            A=randn(n,m);
            A=normc(A);
            
            supp_ind=randperm(m,supp_num);
            supp=zeros(m,1);
            supp(supp_ind)=1;
            x0=zeros(m,1);
            switch dist
                case 1
                    tp=rand>0.5;
                    x0(supp_ind)=tp*(rand(supp_num,1)-2)+(1-tp)*(rand(supp_num,1)+2);
                case 2
                    x0(supp_ind)=randn(supp_num,1);
                case 3
                    x0(supp_ind)=(rand(supp_num,1)>0.5)*2-1;
                case 4
                    x0(supp_ind)=trnd(1,supp_num,1);

            end
            
            b=A*x0;
            
            %% test various algo
            tic;
            [ x, S, r ]=MP(A,b,options);
            err1_MP(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_MP(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(1)=runTimes(1)+toc;
            %
            tic
            [ x, S, r ]=OMP(A,b,options);
            err1_OMP(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_OMP(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(2)=runTimes(2)+toc;
            %
            tic
            [ x, S, r ]=LSOMP(A,b,options);
            err1_LSOMP(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_LSOMP(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(3)=runTimes(3)+toc;
            %
            tic
            [ x, S, r ]=WMP(A,b,options);
            err1_WMP(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_WMP(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(4)=runTimes(4)+toc;
            %
            tic
            [ x, S, r ]=ThMP(A,b,options);
            err1_ThMP(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_ThMP(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(5)=runTimes(5)+toc;
        end
        
    end
    %% plot
    err1_MP_mean=mean(err1_MP,2);
    err2_MP_mean=mean(err2_MP,2);
    
    err1_OMP_mean=mean(err1_OMP,2);
    err2_OMP_mean=mean(err2_OMP,2);
    
    err1_LSOMP_mean=mean(err1_LSOMP,2);
    err2_LSOMP_mean=mean(err2_LSOMP,2);
    
    err1_WMP_mean=mean(err1_WMP,2);
    err2_WMP_mean=mean(err2_WMP,2);
    
    err1_ThMP_mean=mean(err1_ThMP,2);
    err2_ThMP_mean=mean(err2_ThMP,2);
    
    axis_supp_num=1:maxSupp;
    figure;
    %plot(axis_supp_num,err1_MP_mean,axis_supp_num,err1_OMP_mean,axis_supp_num,err1_LSOMP_mean,axis_supp_num,err1_WMP_mean,axis_supp_num,err1_ThMP_mean,'LineWidth',2);
    plot(axis_supp_num,err1_MP_mean,axis_supp_num,err1_OMP_mean,axis_supp_num,err1_LSOMP_mean,axis_supp_num,err1_WMP_mean,'LineWidth',2);
    title('l2 recovery error v.s. Cardinality','FontSize',15,'FontWeight','Bold')
    xlabel('Cardinality of the true solution','FontSize',15);
    ylabel('Average relative l_2 error','FontSize',15);
    legend('MP','OMP','LS-OMP','Weak-MP(t=0.5)','Thresholding','Location','northwest')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/' char(dist_name(dist)) '_l2.eps'] ,'epsc');
    
    figure;
    plot(axis_supp_num,err2_MP_mean,axis_supp_num,err2_OMP_mean,axis_supp_num,err2_LSOMP_mean,axis_supp_num,err2_WMP_mean,axis_supp_num,err2_ThMP_mean,'LineWidth',2);
    title('support error v.s. Cardinality','FontSize',15,'FontWeight','Bold')
    xlabel('Cardinality of the true solution','FontSize',15);
    ylabel('Probability of error in Support','FontSize',15);
    legend('MP','OMP','LS-OMP','Weak-MP(t=0.5)','Thresholding','Location','northwest')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/' char(dist_name(dist)) '_card.eps'] ,'epsc');
    
end
%%
figure;
runTimes_mean=runTimes/maxIter/maxSupp/4;
bar(runTimes_mean);

title('runtimes of the algorithms','FontSize',15,'FontWeight','Bold')
xlabel('Algorithms to compare','FontSize',15);
ylabel('runtime(s)','FontSize',15);
set(gca,'XTickLabel',{'MP','OMP','LS-OMP','Weak-MP(t=0.5)','Thresholding'})
set(gca,'FontSize',10,'FontWeight','Bold');
saveas(gca, ['./eps/runtime.eps' ] ,'epsc');
