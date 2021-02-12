
# identification, C observed, added Lambda matrix, cyclic

# with constraint

G <- matrix(c(1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,1,1),ncol = 5)
G <- G*(abs(rnorm(25)))
G <- G/4

mu <- abs(rnorm(5))

R <- solve(diag(5) - G)
lambda <- R%*%mu
Lambda <- diag(5)
diag(Lambda) <- lambda

C <- R%*%Lambda%*%t(R)

# without constraint

G <- matrix(c(1,1,0,1,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,1,1),ncol = 5)
G <- G*abs(rnorm(25))
G <- G/4

mu <- abs(rnorm(5))

R <- solve(diag(5) - G)
lambda <- R%*%mu
Lambda <- diag(5)
diag(Lambda) <- lambda

C <- R%*%Lambda%*%t(R)

###

R11 <- sqrt(C[1,1])
R11
R[1,1]*sqrt(lambda[1])

R21 <- C[1,2]/(R11)
R21
R[2,1]*sqrt(lambda[1])

R31 <- C[1,3]/(R11)
R31
R[3,1]*sqrt(lambda[1])

R41 <- C[1,4]/(R11)
R41
R[4,1]*sqrt(lambda[1])

gg32 <- R31/R21  
gg32
G[3,2]*(1-G[3,3])^(-1)


lambda[3]*R[3,3]*R[4,3] + G[3,2]*(1-G[3,3])^(-1)*(C[4,2] - lambda[3]*R[4,3]*R[2,3])
lambda[3]*R[3,3]*G[4,3]*(1-G[4,4])^(-1)*R[3,3] + 
  G[3,2]*(1-G[3,3])^(-1)*(C[4,2] - lambda[3]*G[4,3]*(1-G[4,4])^(-1)*R[3,3]*R[2,3])
C[3,4]

lambda[3]*R[3,3]^2 + G[3,2]*(1-G[3,3])^(-1)*(C[3,2] - lambda[3]*R[3,3]*R[2,3])
C[3,3]

(C[3,2] - (C[3,3] - lambda[3]*R[3,3]^2)/(G[3,2]*(1-G[3,3])^(-1)))/(lambda[3]*R[3,3])
R[2,3]

(C[4,2] - (C[3,4] - lambda[3]*R[3,3]^2*G[4,3]*(1-G[4,4]^(-1)))/(G[3,2]*(1-G[3,3])^(-1)))/(lambda[3]*R[3,3]*G[4,3]*(1-G[4,4])^(-1))


gg43<- G[4,3]*(1-G[4,4])^(-1)
(C[4,2] - lambda[3]*C[3,3] + lambda[3]*gg32*C[3,2] - gg32*C[4,2])/(gg32*lambda[3]*(1-gg43))

(C[3,3] - gg32*(C[3,2] - lambda[3]*R[3,3]*R[2,3]))/lambda[3]
R[3,3]^2

(C[3,2] - G[3,2]*(1-G[3,3])^(-1)*(C[2,2] - lambda[3]*R[2,3]^2))/(lambda[3]*R[2,3])
R[3,3]


lambda[3]*((C[3,2] - G[3,2]*(1-G[3,3])^(-1)*(C[2,2] - lambda[3]*R[2,3]^2))/(lambda[3]*R[2,3]))^2 + gg32^2*(C[2,2] - lambda[3]*R[2,3]^2)
C[3,3]

ff <- function(x) 
  C[3,3] - (lambda[3]*((C[3,2] - gg32*(C[2,2] - lambda[3]*x^2))/(lambda[3]*x))^2 + gg32^2*(C[2,2] - lambda[3]*x^2))

plot(ff, 0, 10)

R23 <- uniroot(ff, c(0,1000), tol = .Machine$double.eps/1000)$root
ff(R[2,3])
R23
R[2,3]

R33 <- (C[3,2] - gg32*(C[2,2] - lambda[3]*R23^2))/(lambda[3]*R23)
R33
R[3,3]

(C[3,4] - gg32*C[4,2])/(lambda[3]*R33*(R33-gg32*R23))
R41/R31


###

library(pracma)

O <- diag(5)
O[2:4, 2:4] <- randortho(3)



####################
####################
####################
# identification, C observed, added Lambda matrix, cyclic, big cycle

# with constraint

G <- matrix(c(1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,1,0,0,1,0,0,1,0,1,1),ncol = 5)
G <- G*(abs(rnorm(25)))
G <- G/4

mu <- abs(rnorm(5))

R <- solve(diag(5) - G)
lambda <- R%*%mu
Lambda <- diag(5)
diag(Lambda) <- lambda

C <- R%*%Lambda%*%t(R)

# without constraint

G <- matrix(c(1,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,0,0,1,0,0,1,0,1,1),ncol = 5)
G <- G*abs(rnorm(25))
G <- G/4

mu <- abs(rnorm(5))

R <- solve(diag(5) - G)
lambda <- R%*%mu
Lambda <- diag(5)
diag(Lambda) <- lambda

C <- R%*%Lambda%*%t(R)

###

R11 <- sqrt(C[1,1])
R11
R[1,1]*sqrt(lambda[1])

R21 <- C[1,2]/(R11)
R21
R[2,1]*sqrt(lambda[1])

R31 <- C[1,3]/(R11)
R31
R[3,1]*sqrt(lambda[1])

R41 <- C[1,4]/(R11)
R41
R[4,1]*sqrt(lambda[1])

gg32 <- R31/R21  
gg32
G[3,2]*(1-G[3,3])^(-1)


lambda[3]*R[3,3]*R[4,3] + G[3,2]*(1-G[3,3])^(-1)*(C[4,2] - lambda[3]*R[4,3]*R[2,3])
lambda[3]*R[3,3]*G[4,3]*(1-G[4,4])^(-1)*R[3,3] + 
  G[3,2]*(1-G[3,3])^(-1)*(C[4,2] - lambda[3]*G[4,3]*(1-G[4,4])^(-1)*R[3,3]*R[2,3])
C[3,4]

lambda[3]*R[3,3]^2 + G[3,2]*(1-G[3,3])^(-1)*(C[3,2] - lambda[3]*R[3,3]*R[2,3])
C[3,3]

(C[3,2] - (C[3,3] - lambda[3]*R[3,3]^2)/(G[3,2]*(1-G[3,3])^(-1)))/(lambda[3]*R[3,3])
R[2,3]

(C[4,2] - (C[3,4] - lambda[3]*R[3,3]^2*G[4,3]*(1-G[4,4]^(-1)))/(G[3,2]*(1-G[3,3])^(-1)))/(lambda[3]*R[3,3]*G[4,3]*(1-G[4,4])^(-1))


gg43<- G[4,3]*(1-G[4,4])^(-1)
(C[4,2] - lambda[3]*C[3,3] + lambda[3]*gg32*C[3,2] - gg32*C[4,2])/(gg32*lambda[3]*(1-gg43))

(C[3,3] - gg32*(C[3,2] - lambda[3]*R[3,3]*R[2,3]))/lambda[3]
R[3,3]^2

(C[3,2] - G[3,2]*(1-G[3,3])^(-1)*(C[2,2] - lambda[3]*R[2,3]^2))/(lambda[3]*R[2,3])
R[3,3]


lambda[3]*((C[3,2] - G[3,2]*(1-G[3,3])^(-1)*(C[2,2] - lambda[3]*R[2,3]^2))/(lambda[3]*R[2,3]))^2 + gg32^2*(C[2,2] - lambda[3]*R[2,3]^2)
C[3,3]

ff <- function(x) 
  C[3,3] - (lambda[3]*((C[3,2] - gg32*(C[2,2] - lambda[3]*x^2))/(lambda[3]*x))^2 + gg32^2*(C[2,2] - lambda[3]*x^2))

plot(ff, 0, 10)

R23 <- uniroot(ff, c(0,1000), tol = .Machine$double.eps/1000)$root
ff(R[2,3])
R23
R[2,3]

R33 <- (C[3,2] - gg32*(C[2,2] - lambda[3]*R23^2))/(lambda[3]*R23)
R33
R[3,3]

id1 <- (C[3,4] - gg32*C[4,2])/(lambda[3]*R33*(R33-gg32*R23))


###
# id R41/R31

# id R41

gg14 <- G[1,4]/(1-G[1,1])

ff <- function(x) 
  C[3,3] - (lambda[3]*((C[3,2] - gg32*(C[2,2] - lambda[3]*x^2))/(lambda[3]*x))^2 + gg32^2*(C[2,2] - lambda[3]*x^2))

ff11 <- function(x) 
  C[1,1] - (lambda[1]*((C[1,4] - gg14*(C[4,4] - lambda[1]*x^2))/(lambda[1]*x))^2 + gg14^2*(C[4,4] - lambda[1]*x^2))


plot(ff11, 0, 10)

R41 <- uniroot(ff11, c(0,1000), tol = .Machine$double.eps/1000)$root
ff11(R[4,1])
R41
R[4,1]

# id R11

R11 <- (C[1,4] - gg14*(C[4,4] - lambda[1]*R41^2))/(lambda[1]*R41)
R11
R[1,1]

# id R31


ff31 <- function(x){
  C[3,1] - (lambda[1]*x*R11 + gg14*(C[3,4] - lambda[1]*x*R41))}


plot(ff31, 0, 10)

R31 <- uniroot(ff31, c(0,1000), tol = .Machine$double.eps/1000)$root
ff31(R[3,1])
R31
R[3,1]

id2 <- R41/R31
R[4,1]/R[3,1]


################

library(SEMID)

Gtmp <- G
diag(Gtmp) <- 0
Otmp <- diag(nrow(G))

Ggraph <- MixedGraph(L = t(Gtmp > 0), O = matrix(0, nrow = nrow(Gtmp), ncol = ncol(Gtmp)))
rm(Gtmp)
plot(Ggraph)

ed <- edgewiseID(Ggraph)
t(ed$identifier((R)%*%t(R))$Lambda)*(1-G[2,2])
G
tmp <- diag(5)
diag(tmp) <- 1 - diag(G)
tmp%*%t(ed$identifier((R)%*%t(R))$Lambda) - G

ed$identifier((R)%*%t(R))$Omega
1/tmp^2
rm(tmp)

# check first that this is identifiable, actually

idG <- function(C, g){
  
  G0 <- matrix(0, nrow = nrow(C), ncol = ncol(C))
  
  ed <- edgewiseID(g)
  
  diag(G0) <- sqrt(1/diag(ed$identifier(t(C))$Omega))
  
  tmp <- G0%*%t(ed$identifier(t(C[-5,-5]))$Lambda)
  diag(tmp) <- 1 - diag(G0)
  
  tmp
}



Gtmp <- G
diag(Gtmp) <- 0
Otmp <- matrix(0, nrow = nrow(Gtmp), ncol = ncol(Gtmp))
Otmp[2,4] <- Otmp[4,2] <- 1
Otmp <- Otmp[-5,-5]

Ggraph <- MixedGraph(L = t(Gtmp[-5,-5] > 0), O = Otmp)
rm(Gtmp)
plot(Ggraph)

idG(C[-5,-5], Ggraph)
G[-5,-5]

G2 <- G
G2[1,1] <- G2[1,1]*2
G2[2,1] <- G2[2,1]/2
G2[2,2] <- G2[2,2]*2

G2[1,1] <- G2[1,1] - .1

ss <- matrix(0, nrow = 5, ncol = 5)
ss[1,1] <- .1

R%*%ss

R2 <- solve(diag(5) - G2)

R2 <- R
R2[4,5]

C2 <- R2%*%t(R2)
C2
C

G <- matrix(0, nrow= 2,ncol=2)
G[1,1] <- .1
G[2,1] <- .3
G[2,2] <- .4
R <- solve(diag(2)-G)
R%*%t(R)

G0 <- matrix(0, nrow = nrow(C[-5,-5]), ncol = ncol(C[-5,-5]))

ed <- edgewiseID(g)

1-1/sqrt((ed$identifier(t(C[-5,-5]))$Omega))

Om <- ed$identifier(t(C[-5,-5]))$Omega
e <- eigen(Om)
e$vectors%*%diag(e$values)%*%t(e$vectors)

solve(e$vectors%*%sqrt(diag(1/e$values)))

1-1/sqrt(diag(Om))

diag(G0) <- 1-diag(G)[-5]

diag(G)

tmp <- G0%*%t(ed$identifier(t(C[-5,-5]))$Lambda)
diag(tmp) <- 1 - diag(G0)

G0%*%t(ed$identifier(t(C[-5,-5]))$Lambda) + diag(4) - G0
G[-5,-5]

t(ed$identifier(t(C[-5,-5]))$Lambda)
solve(G0)%*%(G[-5,-5] - diag(4) + G0)
solve(G0)%*%G[-5,-5] - solve(G0) + diag(4)


tmp

rr <- solve(diag(4) - G0)%*%tt%*%solve(diag(4) - G0)
diag(solve(diag(4) - G0))^2

sqrt(diag(rr))

# task: id diag of G using Omega




# not done: Lambda (mean from Hawkes process should be added to equations)

###
# small example

G <- matrix(abs(rnorm(4)), ncol = 2)
G[2,1] <- 0
R <- solve(diag(2)-  G)
C <- t(R)%*%R


Gtmp <- G
diag(Gtmp) <- 0
Otmp <- diag(nrow(G))

Ggraph <- MixedGraph(L = t(Gtmp > 0), O = matrix(0, nrow = nrow(Gtmp), ncol = ncol(Gtmp)))
rm(Gtmp)
plot(Ggraph)

ed <- edgewiseID(Ggraph)
t(ed$identifier((G)%*%t(G))$Lambda)


# make fake nodes

G1 <- matrix(0, ncol = 4, nrow = 4)

G1[1,2] <- G[1,1]
G1[2,1] <- G[1,1]
G1[3,4] <- G[2,2]
G1[4,3] <- G[2,2]
G1[1,3] <- G1[1,4] <- G1[2,3] <- G1[2,4] <- G[1,2]
G1[3,1] <- G1[3,2] <- G1[4,1] <- G1[4,2] <- G[2,1]


Gtmp <- G1
diag(Gtmp) <- 0
Otmp <- diag(nrow(G))

Ggraph <- MixedGraph(L = t(Gtmp > 0), O = matrix(0, nrow = nrow(Gtmp), ncol = ncol(Gtmp)))
rm(Gtmp)
plot(Ggraph)

ed <- edgewiseID(Ggraph)
t(ed$identifier(t(G1)%*%G1)$Lambda)

