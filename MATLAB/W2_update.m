function [X3H2t,W2H2H2t] = W2_update(Z,W2,H2)
    
    H2t = H2';   
    X3H2t = Z*H2t;
    H2H2t = H2*H2t;
    W2H2H2t = W2*H2H2t;
  
end