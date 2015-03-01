function [ x, S, r ]=OMP( A,b,options )
%
% input: A normalized!
% modified version, recursive LLS to speed up

min_err=options.min_error;
b0=b;
[n,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);
Q=zeros(n,m);
iter=1;
while norm(b)^2 > min_err
    
    max_corr=0;
    min_i=0;
    min_q=0;
    
    I=1:m;
    
    for i=I(S==0)
        
        if(iter>1)
            tq=projectQperp(Q(:,1:iter-1),A(:,i));
        else
            tq=A(:,i); % A is assumed normalized
        end
        tcorr=abs(tq'*b);
        if (tcorr>max_corr)
            min_i=i;
            max_corr=tcorr;
            min_q=tq;
        end
    end
    Q(:,iter)=min_q;
    b=b-min_q*min_q'*b;
    S(min_i)=1;
    iter=iter+1;
end
r=norm(b)/norm(b0);
x(S==1)=A(:,S==1)\b0;
end

