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
    coefs=A'*b;
    [~,ind]=max(abs(coefs)); % assuming only one max exists
    S(ind)=1;
    if (iter>1)
        q=projectQperp(Q(:,1:iter-1),A(:,ind));
        
    else
        q=A(:,ind); % A is assumed normalized
    end
    b=b-q*q'*b;
    Q(:,iter)=q;
    iter=iter+1;
end

x(S==1)=A(:,S==1)\b0;
r=norm(b)/norm(b0);
end

