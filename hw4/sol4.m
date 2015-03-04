clc;
close all;
%% generating A and normalize
n=30;m=50;

%% eval random x for uniform, gaussian , {+1,-1}, cauchy
maxIter=1000;
maxSupp=10;

options.min_error=1e-4;
options.t=0.5;
runTimes=zeros(1,7);
    
    
    err1_MP=zeros(maxSupp,maxIter);
    err2_MP=zeros(maxSupp,maxIter);
    
    err1_OMP=zeros(maxSupp,maxIter);
    err2_OMP=zeros(maxSupp,maxIter);
    
    err1_LSOMP=zeros(maxSupp,maxIter);
    err2_LSOMP=zeros(maxSupp,maxIter);
    
    err1_WMP=zeros(maxSupp,maxIter);
    err2_WMP=zeros(maxSupp,maxIter);
    
    err1_lasso=zeros(maxSupp,maxIter);
    err2_lasso=zeros(maxSupp,maxIter);
    
    err1_rel1=zeros(maxSupp,maxIter);
    err2_rel1=zeros(maxSupp,maxIter);
    
    err1_rel2=zeros(maxSupp,maxIter);
    err2_rel2=zeros(maxSupp,maxIter);
    
    for supp_num=1:maxSupp
        disp(['-------- supp= ' num2str(supp_num) ' / ' num2str(maxSupp)])
        for iter=1:maxIter
            if (mod(iter,100)==0)
                    disp(['------------ iter= ' num2str(iter) ' / ' num2str(maxIter)])
            end

            A=randn(n,m);
            A=normc(A);
            
            supp_ind=randperm(m,supp_num);
            supp=zeros(m,1);
            supp(supp_ind)=1;
            x0=zeros(m,1);
                    tp=rand>0.5;
                    x0(supp_ind)=tp*(rand(supp_num,1)-2)+(1-tp)*(rand(supp_num,1)+2);

            
            
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
            cvx_begin quiet
                variable x(m)
                minimize(norm(x,1));
                subject to
                     A*x == b;
            cvx_end
            
            %X=lasso(A,b,'Lambda',1e-10);

            %x=X(:,1);

            S=(abs(x)>1e-4);
            err1_lasso(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_lasso(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(5)=runTimes(5)+toc;
            
            %
            tic
            x=lp_re(A,b,1);
            S=(abs(x)>1e-4);
            err1_rel1(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_rel1(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(6)=runTimes(6)+toc;
            %
            tic
            x=lp_re(A,b,2);
            S=(abs(x)>1e-4);
            err1_rel2(supp_num,iter)=norm(x-x0)^2/norm(x0)^2;
            err2_rel2(supp_num,iter)=(max(sum(S),supp_num)-sum(S.*supp))/max(sum(S),supp_num);
            runTimes(7)=runTimes(7)+toc;
            
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
    
    err1_lasso_mean=mean(err1_lasso,2);
    err2_lasso_mean=mean(err2_lasso,2);
    
    err1_rel1_mean=mean(err1_rel1,2);
    err2_rel1_mean=mean(err2_rel1,2);
    
    err1_rel2_mean=mean(err1_rel2,2);
    err2_rel2_mean=mean(err2_rel2,2);
    
    axis_supp_num=1:maxSupp;
    figure;
    %plot(axis_supp_num,err1_MP_mean,axis_supp_num,err1_OMP_mean,axis_supp_num,err1_LSOMP_mean,axis_supp_num,err1_WMP_mean,axis_supp_num,err1_ThMP_mean,'LineWidth',2);
    %plot(axis_supp_num,err1_MP_mean,'-',axis_supp_num,err1_OMP_mean,'--',axis_supp_num,err1_LSOMP_mean,':',axis_supp_num,err1_WMP_mean,'-.','LineWidth',2);
    plot(axis_supp_num,err1_MP_mean,'-',axis_supp_num,err1_OMP_mean,'--',axis_supp_num,err1_LSOMP_mean,':',axis_supp_num,err1_WMP_mean,'-.',axis_supp_num,err1_lasso_mean,'+-',axis_supp_num,err1_rel1_mean,'o',axis_supp_num,err1_rel2_mean,'^','LineWidth',2);
    
    title('l2 recovery error v.s. Cardinality','FontSize',15,'FontWeight','Bold')
    xlabel('Cardinality of the true solution','FontSize',15);
    ylabel('Average relative l_2 error','FontSize',15);
    legend('MP','OMP','LS-OMP','Weak-MP(t=0.5)','lasso','re-l1','re-l2','Location','northwest')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/' 'l2.eps'] ,'epsc');
    
    figure;
    plot(axis_supp_num,err2_MP_mean,'-',axis_supp_num,err2_OMP_mean,'--',axis_supp_num,err2_LSOMP_mean,':',axis_supp_num,err2_WMP_mean,'-.',axis_supp_num,err2_lasso_mean,'+-',axis_supp_num,err2_rel1_mean,'o',axis_supp_num,err2_rel2_mean,'^','LineWidth',2);
    title('support error v.s. Cardinality','FontSize',15,'FontWeight','Bold')
    xlabel('Cardinality of the true solution','FontSize',15);
    ylabel('Probability of error in Support','FontSize',15);
    legend('MP','OMP','LS-OMP','Weak-MP(t=0.5)','lasso','re-l1','re-l2','Location','northwest')
    set(gca,'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['./eps/'  'card.eps'] ,'epsc');
    

%%
figure;
runTimes_mean=runTimes/maxIter/maxSupp;
bar(runTimes_mean);

title('runtimes of the algorithms','FontSize',15,'FontWeight','Bold')
xlabel('Algorithms to compare','FontSize',15);
ylabel('runtime(s)','FontSize',15);
set(gca,'XTickLabel',{'MP','OMP','LS-OMP','Weak-MP(t=0.5)','lasso','re-l1','re-l2'})
set(gca,'FontSize',10,'FontWeight','Bold');
saveas(gca, ['./eps/runtime.eps' ] ,'epsc');
