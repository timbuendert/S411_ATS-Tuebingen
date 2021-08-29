function [D,A,B_0,C] = Cholesky_decomposition(vcv)

  
  pp = chol(vcv);
  P = pp';                       % Lower triangular matrix
  C = diag(P).*eye(size(vcv,2));
  D = C*C';
  A = P/C;
  B_0 = inv(A);
  

end