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
    
    min_b=inf;
    min_i=0;
    min_x=0;
    
    I=1:m;
    
    for i=I(S==0)
        tS=S;
        tS(i)=1;
        tA=A(:,tS==1);
        tx=tA\b0;
        if (norm(b0-tA*tx)<norm(min_b))
            min_i=i;
            min_b=b0-tA*tx;
            min_x=tx;
        end
    end
    
    b=min_b;
    S(min_i)=1;
    x(S==1)=min_x;
end
r=norm(b)/norm(b0);

end

