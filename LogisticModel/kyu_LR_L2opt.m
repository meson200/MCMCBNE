function [w,conv] = kyu_LR_L2opt(X,Y,lambda,maxiter)

% optimize weights of a L2 logistic function 
% at a fixed lambda (L2 term coeff.)
% that minimizes least-squares
% using iterative Newton-Raphson method
% (Pattern recognition (Bishop, 2006) chapter 4.3)
%
% input
% X: data matrix (without the intercept!)
% Y: class column
% lambda: L2 regularization coeff.
% maxiter: maximum iteration 
%
% output
% w: coefficients for variables (including an intercept)
% conv: if convergence was reached (1) or not (0)

% add the intercept
X = [X ones(size(X,1),1)];

m = size(X,1);
t = 1;
wdel = 1000;
errt = [];
npar = size(X,2);
w = 0.1*ones(npar,1);
while wdel > 0.001 && t<maxiter
        t = t+1;
        h_trntrn = h_log(w,X);
        err = Y - h_trntrn;
        sumerr = sum(abs(err));

        %make R matrix
        R = eye(m);
        for kk = 1:m
            R(kk,kk) = h_trntrn(kk)*(1-h_trntrn(kk));
        end
        % Newton-Raphson update formula with L2 regularization 
        wnew = w + (X'*R*X + lambda*eye(npar))\(X'*err - lambda*w);
        %wnew = w + alpha*L_del(err,X_train);
        wdel = norm(wnew-w,2);
        errt = [errt sumerr];
        w = wnew;
end
if t<maxiter 
   conv = 1;
else
   conv = 0;
end

function hval = h_log(w,X)
% logit link function
l = size(X,1);
hval = zeros(l,1);
for k = 1:l
    hval(k) = 1/(1+exp(-X(k,:)*w));
end