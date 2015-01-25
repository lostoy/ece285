function [ x, S, r ]=ThMP( A,b,options )
%
% input: A normalized!
% modified version, recursive LLS to speed up
min_err=options.min_error;


b0=b;

[n,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);

coefs=A'*b;
[~,inds]=sort(abs(coefs),'descend');
Q=zeros(n,m);
iter=1;

while norm(b)^2 > min_err
    if (iter>1)
       q=projectQperp(Q(:,1:iter-1),A(:,inds(iter))); 
    else
        q=A(:,inds(1));
    end
    
    b=b-q*q'*b;
    Q(:,iter)=q;
    S(inds(iter))=1;
    iter=iter+1;
end


r=norm(b)/norm(b0);
x(S==1)=A(:,S==1)\b0;
end

