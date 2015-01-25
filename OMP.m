function [ x, S, r ]=OMP( A,b,options )
%
% input: A normalized!

min_err=options.min_error;
b0=b;
[n,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);
%Q=zeros(n,m); % TODO: implicit Q of QR

while norm(b)^2 > min_err
    coefs=A'*b;
    [~,ind]=max(abs(coefs)); % assuming only one max exists
    S(ind)=1;
    
    x1=A(:,S==1)\b0;
    b=b0-A(:,S==1)*x1; %TODO: implicit Q of QR
    x(S==1)=x1;
    
end
r=norm(b)/norm(b0);


end

