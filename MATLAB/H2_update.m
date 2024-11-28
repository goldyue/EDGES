function [W1tW1H2,W1tX2,W2tX3,W2tW2H2,ekkH2] = H2_update(W1,W2,H2,Z,Y,K)
    
    W1t = W1';
    W2t = W2'; 
    %H2L2 = H2*L2_final;
    W1tX2 = W1t*Y;
    W2tX3 = W2t*Z; 
    W1tW1H2 = W1t*W1*H2;
    W2tW2H2 = W2t*W2*H2; 
    K=20;
    ekk = eye(K,K);
    ekkH2 = ekk*H2;
    
end