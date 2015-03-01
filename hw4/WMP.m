function [ x, S, r ]=WMP( A,b,options )
%
% input: A normalized!

min_err=options.min_error;
t=options.t;
b0=b;
[~,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);

while norm(b)^2 > min_err
    coefs=A'*b;
    ind=find(abs(coefs)>t*norm(b),1,'first');
    if(isempty(ind))
        [~,ind]=max(abs(coefs));
    end
    S(ind)=1;
    coef=coefs(ind);
    b=b-coef*A(:,ind);
    x(ind)=x(ind)+coef;
end
r=norm(b)/norm(b0);

end

