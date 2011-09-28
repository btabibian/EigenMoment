function [A,B]=construct_AB(m,c)
A=zeros(m,m);
for i=0:m-1
    for j=0:m-1
        syms x1 x2
        f=x1^j*x2^i*exp(-1*c*(x1-x2)^2);
        g=inline(vectorize(f),'x1','x2');
        A(i+1,j+1)= dblquad(g,-1,1,-1,1);
    end
end
B=zeros(m,m);
for i=0:m-1
    for j=0:m-1
        B(i+1,j+1)= 1/(i+j+1)*(1^(i+j+1)-(-1)^(i+j+1));
    end
end



