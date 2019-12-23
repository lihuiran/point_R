# findpoints
#
# This is an example function named 'findpoints'
# which return label.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

label <- function(X,k,branch,tip) {
  m <- nrow(X)
  n <- ncol(X)
  label <- matrix(0,nrow=n,ncol=1)
  X = t(X)
  dis <- dis(X)

  D <- matrix(0,nrow=n,ncol=n)
  I <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    D[,i] <- sort(dis[,i])
    I[,i] <- order(dis[,i])
  }

  X <- t(X)
  sigma <- matrix(0,nrow=m,ncol=n)
  d <- matrix(0,nrow=1,ncol=n)
  for(i in 1:n){
    Ni <- X[,I[1:k,i]]
    centeri <- rowMeans(Ni)
    CNi <- Ni - matrix(rep(centeri,k),nrow=m,ncol=k)
    S <- svd(CNi)
    sigma[,i] <- S$d
    d[i] <- sigma[1,i]/sigma[2,i]
  }

  I3 <- find(d<branch)
  I12 <- find(d>=branch)

  vector <- array(0,dim=c(m,k,length(I12)))
  for(i in 1:length(I12)){
    Ni <- X[,I[1:k,I12[i]]]
    for (j in 2:k){
      vector[,j,i] <- Ni[,j]-X[,I12[i]]
    }
  }

  cosine <- array(0,dim=c(k,k,length(I12)))
  for(i in 1:length(I12)){
    for (p in 2:k){
      for (q in 2:k){
        cosine[,,i][p,q] <- t(vector[,p,i])%*%vector[,q,i]/sqrt(sum((vector[,q,i])^2))/sqrt(sum((vector[,p,i])^2))
      }
    }
  }


  negative <- matrix(0,nrow=1,ncol=length(I12))
  positive <- matrix(0,nrow=1,ncol=length(I12))
  rate <- matrix(0,nrow=1,ncol=length(I12))
  for(i in 1:length(I12)){
    cosi <- cosine[,,i]
    negative[i] <- length(find(cosi<0))/2
    positive[i] <- (k-1)*(k-2)-negative[i]
    rate[i] <- negative[i]/((k-1)*(k-2)/2)
  }

  I1 <- find(rate<tip)
  I1 <- I12[I1]
  label[I3] <- 3
  label[I1] <- 1
  I2 <- find(label==0)
  label[I2] <- 2

  return(label)
}
