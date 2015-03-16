function x=SBL(A,b)
[n,m]=size(A);
gammas=eye(m);
sigma=1;

maxIter=400;
for iter=1:maxIter
    mu=gammas*A'*(A*gammas*A'+sigma*eye(n))^-1*b;
    P=gammas-gammas*A'*(A*gammas*A'+sigma*eye(n))^-1*A*gammas;
    
    gammas0=gammas;
    sigma=(b'*b-2*mu'*A'*b+trace(A'*A*P))/n;
    gammas=diag(diag(mu*mu'+P));
    if (norm(gammas-gammas0,'fro')<1e-2*norm(gammas0,'fro'))
        break;
    end
end

x=mu;
end