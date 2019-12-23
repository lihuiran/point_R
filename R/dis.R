
dis <- function(X){
m <- nrow(X)
n <- ncol(X)
A <- X%*%t(X)
l <- t(diag(A))
l <- t(l)
Dis <- sqrt(matrix(rep(l,m),nrow=m,ncol=m)+matrix(rep(t(l),each=m),nrow=m,ncol=m)-2*A)
}
