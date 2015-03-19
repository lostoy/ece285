function x=reweightedL1(A,b,x)
[n,m]=size(A);

maxIter=5000;
epsilon=1e-6;

for iter=1:maxIter
    W=diag(1./(abs(x)+epsilon));
    cvx_begin quiet
        variable x1(m)
        minimize(norm(W*x1,1));
        subject to
        A*x1 == b;
    cvx_end
    x=x1;
    if (norm(x-x1)<1e-6*norm(x))
        break;
    end
end

end