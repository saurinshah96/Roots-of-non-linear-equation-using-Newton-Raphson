clc
close all

syms x1 x2 x3
f1(x1,x2,x3) = x1^2 - 2*x1 + x2^2 - x3 + 1;
f2(x1,x2,x3) = x1*x2*x2 - x1 - 3*x2 + x2*x3 + 2;
f3(x1,x2,x3) = x1*x3*x3 - 3*x3 + x2*x3*x3 + x1*x2;

J(x1,x2,x3) = jacobian([f1(x1,x2,x3),f2(x1,x2,x3), f3(x1,x2,x3)]);
%J_det = det(J(0,0,0));

 inv_J = inv(J);

q1 = 1;
q2 = 2;
q3 = 3;
  
for i = 0:0
    r0 = [q1;q2;q3];
    Jac = inv_J(q1,q2,q3);
    v1 = f1(q1,q2,q3);
    v2 = f2(q1,q2,q3);
    v3 = f3(q1,q2,q3);
    Fx = [v1; v2; v3];
    Inner = Jac*Fx;
    Ans = r0 - Inner;
    Ans = vpa(Ans)
    q1 = double(Ans(1,1));
    q2 = double(Ans(2,1));
    q3 = double(Ans(3,1));

end



