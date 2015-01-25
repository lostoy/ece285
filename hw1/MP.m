function [ x, S, r ]=MP( A,b,options )
%
% input: A normalized!

min_err=options.min_error;
b0=b;
[~,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);

while norm(b)^2 > min_err
    coefs=A'*b;
    [~,ind]=max(abs(coefs)); % assuming only one max exists
    coef=coefs(ind);
    S(ind)=1;
    b=b-coef*A(:,ind);
    x(ind)=x(ind)+coef;
    
end
r=norm(b)/norm(b0);

end

