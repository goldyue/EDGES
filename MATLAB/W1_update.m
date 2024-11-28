function [X1H1t,X2H2t,W1H1H1t,W1H2H2t] = W1_update(W1,H1,H2,X,Y)

    H1t = H1';
    H2t = H2';
    X1H1t = X*H1t;
    X2H2t = Y*H2t;
    H2H2t = H2*H2t;
    H1H1t = H1*H1t;
    W1H1H1t = W1*H1H1t;
    W1H2H2t = W1*H2H2t;
  
end