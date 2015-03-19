function x=reweightedL2(A,b,x)
[n,m]=size(A);

maxIter=50000;
epsilon=1e-6;

for iter=1:maxIter
    W=diag(1./((x.^2+epsilon)));
    cvx_begin quiet
        variable x1(m)
        minimize((sum(W*(x1.^2))));
        subject to
        A*x1 == b;
    cvx_end
    x=x1;
    if (norm(x-x1)<1e-7*norm(x))
        break;
    end
end

end