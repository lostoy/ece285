function [ x, S, r ]=ThMP( A,b,options )
%
% input: A normalized!

min_err=options.min_error;


b0=b;

[~,m]=size(A);
S=zeros(m,1);
x=zeros(m,1);
k=1;
coefs=A'*b;
[~,inds]=sort(abs(coefs),'descend');
while norm(b)^2 > min_err
    S=zeros(m,1);
    x=zeros(m,1);
    
    S(inds(1:k))=1;
    x1=A(:,S==1)\b0;
    b=b0-A(:,S==1)*x1;
    x(S==1)=x1;
    
    k=k+1;
end


r=norm(b)/norm(b0);

end

