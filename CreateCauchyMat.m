function mat = CreateCauchyMat(n, m)
mat = sparse(n+m-1+n,n+m-1+n);
for x=1:n+m-1+n
    mat(x, x:(x+n-1)) = [n:-1:1];
end
mat = mat(1:(n+m-1), n:(n+m-1));


