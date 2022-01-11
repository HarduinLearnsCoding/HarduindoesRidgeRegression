function [kernelmatrix] = kernelhw05(x,y,N,gamma)
for m=1:N
    for n=1:N
        kernelmatrix(m,n)=exp(-gamma*(x(n)-y(m))^2);
    end
end
end

