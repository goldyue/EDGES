function [W1tX1,W1tW1H1,H1L1,ekkH1] = H1_update(W1,H1,X,L1_final,K)
    
    W1t = W1';
    H1L1 = H1*L1_final;
    W1tX1 = W1t*X;
    W1tW1H1 = W1t*W1*H1;
    ekk = eye(K,K);
    ekkH1 = ekk*H1;
   
    
end