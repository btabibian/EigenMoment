function [W]=EM(A,B,k)
[V_B,D_B]=eig(B);
V_B_selected=V_B(:,1:k);
D_B_selected=D_B(1:k,1:k);
W1=V_B_selected*D_B_selected^(-1/2);
W1'*B*W1;
W1T_A_W1=W1'*A*W1;
[V_W1T_A_W1,D_W1T_A_W1]=eig(W1T_A_W1);
W2=V_W1T_A_W1;
W=W1*W2;
