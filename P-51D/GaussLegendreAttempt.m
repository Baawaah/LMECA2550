%Gauss Legendre Test
% pas n
%interval

A = [0,4];
n = 100;
h = 1/n;
X = [0.7746 0 -0.7746];
W = [0.5556 0.8889 0.5556];

L = A(2)-A(1);
B = linspace(A(1),A(2),n);

acc = 0;
fct = @(x) 2*x;
for i = 1 : length(B)-1
  alpha = (B(i+1)-B(i))/2;
  beta  = (B(i+1)+B(i))/2; 
  acc =  alpha * ( W(1)*fct( alpha*X(1) + beta) + W(2)*fct( alpha*X(2) + beta) + W(3)*fct( alpha*X(3) + beta) ) + acc;
end
acc