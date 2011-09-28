function [EM]= EigenMoments(m,c,k,Moments)
[A,B]=construct_AB(m,c);
[W]= EM_construct(A,B,k);
EM=transform_image_moments(W,Moments);
function [W]=EM_construct(A,B,k)
[V_B,D_B]=eig(B);
V_B_selected=V_B(:,1:k);
D_B_selected=D_B(1:k,1:k);
W1=V_B_selected*D_B_selected^(-1/2);
%W1'*B*W1;
W1T_A_W1=W1'*A*W1;
[V_W1T_A_W1,D_W1T_A_W1]=eig(W1T_A_W1);
W2=V_W1T_A_W1;
W=W1*W2;
function [O]=transform_image_moments(W,m)
O=W'*m*W;
function [A,B]=construct_AB(m,c)
A=zeros(m,m);
for i=0:m-1
    for j=0:m-1
        syms x1 x2
        f=x1^i*x2^j*exp(-1*c*(x1-x2)^2);
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
