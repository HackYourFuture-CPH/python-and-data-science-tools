###
# examples for hawkesEq paper

library(SEMID)

#

rm(list = ls())

####
# Verma-graph with small cycle

G <- matrix(c(1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,1,1),ncol = 5)
G <- G*(abs(rnorm(25)))
G <- G/4
eigen(G)$values

mu <- abs(rnorm(5))

R <- solve(diag(5) - G)
lambda <- R%*%mu
Lambda <- diag(5)
diag(Lambda) <- lambda

C <- R%*%Lambda%*%t(R)

# graph

Gtmp <- G[1:4,1:4]
diag(Gtmp) <- 0
Otmp <- matrix(0, nrow = 4, ncol = 4)
Otmp[2,4] <- Otmp[4,2] <- 1

Ggraph <- MixedGraph(L = 1*t(Gtmp > 0), O = Otmp)
rm(Gtmp, Otmp)
plot(Ggraph)

ed <- edgewiseID(Ggraph)
Gt <- t(ed$identifier(C[1:4,1:4])$Lambda) # Gtilde
Th <- ed$identifier(C[1:4,1:4])$Omega # Theta

# identification of normalized G

tmp <- diag(5)
diag(tmp) <- 1/(1-diag(G))
normG <- tmp%*%G
rm(tmp)
diag(normG) <- 0
max(abs(Gt - normG[1:4,1:4]))

# identification of (unnormalized) G

Gi <- matrix(NA, nrow = 4, ncol = 4)

diag(Gi) <- 1 - sqrt(lambda[-5]/diag(Th))

tmp <- matrix(0, nrow = 4, ncol = 4)
diag(tmp) <- diag(Gi)

Gi <- (diag(4)-tmp)%*%Gt
diag(Gi) <- diag(tmp)

all.equal(Gi[1,1], G[1,1])
all.equal(Gi[3,2], G[3,2])
all.equal(Gi[3,3], G[3,3])


############################
# example to show that G21, G22, G43, G44, G42 are not identifiable from integrated covariance

Th.f <- function(G,Lambda){
  O <- 1:4
  U <- 5
  D <- diag(1/(1-diag(G)))
  Gtilde <- D%*%G - diag(5)
  diag(Gtilde) <- 0
  (D%*%Lambda%*%D)[O,O] + (diag(5) - Gtilde)[O,U,drop=FALSE]%*%solve(diag(5) - Gtilde)[U,U,drop=FALSE]%*%(D%*%Lambda%*%D)[U,U,drop=FALSE]%*%t(solve(diag(5) - Gtilde)[U,U,drop=FALSE])%*%t((diag(5) - Gtilde))[U,O,drop=FALSE]
}


a1 <- .1
b1 <- 4

d2 <- Gtilde[2,5]^2*Lambda[5,5]/(1-G[5,5])^2 - (a1*Gtilde[2,5])^2*b1*Lambda[5,5]/(1-G[5,5])^2
d4 <- Gtilde[4,5]^2*Lambda[5,5]/(1-G[5,5])^2 - ((1/(a1*b1))*Gtilde[4,5])^2*b1*Lambda[5,5]/(1-G[5,5])^2

G2 <- G
G2[2,2] <- 1 - 1/sqrt(1/(1-G[2,2])^2 + d2/Lambda[2,2])
G2[4,4] <- 1 - 1/sqrt(1/(1-G[4,4])^2 + d4/Lambda[4,4])
G2[2,5] <- (1-G2[2,2])*a1*G[2,5]/(1-G[2,2]) #multiplication needs to be in Gtilde!
G2[4,5] <- (1-G2[4,4])*(1/(a1*b1))*G[4,5]/(1-G[4,4])

Lambda2 <- Lambda
Lambda2[5,5] <- b1*Lambda[5,5]

G2[2,1] <- (1-G2[2,2])*G[2,1]/(1-G[2,2])
G2[2,4] <- (1-G2[2,2])*G[2,4]/(1-G[2,2])
G2[4,3] <- (1-G2[4,4])*G[4,3]/(1-G[4,4])

R <- solve(diag(5) - G)
R2 <- solve(diag(5) - G2)

C <- R%*%Lambda%*%t(R)
C2 <- R2%*%Lambda2%*%t(R2)

max(abs(C[1:4,1:4] - C2[1:4,1:4]))

min(G)
min(G2)

min(Lambda)
min(Lambda2)

abs(eigen(G)$values)
abs(eigen(G2)$values)

G
G2

# id'ed
all.equal(G[1,1], G2[1,1])
all.equal(G[3,2], G2[3,2])
all.equal(G[3,3], G2[3,3])

# not id'ed
all.equal(G[2,1], G2[2,1])
all.equal(G[2,2], G2[2,2])
all.equal(G[4,3], G2[4,3])
all.equal(G[4,4], G2[4,4])
all.equal(G[2,4], G2[2,4])

###################
# example to show constraint in small cycle example
# using Theorem 1 and Lemma 1 of Chen et al. (2014)

R[4,1]/R[3,1]
normG[4,3]
C[1,4]/C[1,3]
(C[3,4] - normG[3,2]*C[2,4])/(C[3,3] - normG[3,2]*C[2,3])
((diag(5) - normG)%*%C)[3,4]/((diag(5) - normG)%*%C)[3,3]

###################
###################
# find example to show that nothing is identifiable in big cycle example (and no constraints) without additional info (or just that this specific equality constraint does not hold, the equality holds, but is not identifiable from the integrated covariance)

rm(list = ls())

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

Gtmp <- G[1:4,1:4]
diag(Gtmp) <- 0
Otmp <- matrix(0, nrow = 4, ncol = 4)
Otmp[2,4] <- Otmp[4,2] <- 1

Ggraph <- MixedGraph(L = 1*t(Gtmp > 0), O = Otmp)
rm(Gtmp, Otmp)
plot(Ggraph)

# semID is inconclusive, try different criteria, or find example?
# check Chen's papers?

graphID.nonHtcID(L = 1*t(Gtmp > 0), O = Otmp)
ed <- semID(Ggraph, testGenericNonID = TRUE)
ed <- generalGenericID(Ggraph, idStepFunctions = list(htcIdentifyStep,
                                                ancestralIdentifyStep,
                                                edgewiseIdentifyStep,
                                                trekSeparationIdentifyStep))
Gt <- t(ed$identifier(C[1:4,1:4])$Lambda) # Gtilde
Th <- ed$identifier(C[1:4,1:4])$Omega # Theta

###

D <- diag(1/(1-diag(G)))
Gtilde <- D%*%G - diag(5)
diag(Gtilde) <- 0

x<-.1
OO <- matrix(c(1,x,x*Gtilde[3,2],x*Gtilde[3,2]*Gtilde[4,3],0,
               x*Gtilde[4,3]*Gtilde[4,1],1,x,x*Gtilde[4,3],0,
               x*Gtilde[1,4],0,1,x,0,
               x,x*Gtilde[2,1],0,1,0,
               0,0,0,0,0), 
             nrow = 5)

Gtilde2 <- (diag(5)-Gtilde)%*%OO





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

##########################
##########################
##########################
# id using dynamical info (known parameters)

# identification, C observed, cyclic, big cycle

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

#

tmp <- diag(5)
diag(tmp) <- 1/(1-diag(G))
normG <- tmp%*%G
rm(tmp)
diag(normG) <- 0

D <- diag(5)
diag(D) <- 1-diag(G)

max(abs(C - solve(diag(5) - normG)%*%solve(D)%*%
          Lambda%*%solve(D)%*%t(solve(diag(5) - normG))))

# defining auxiliary variable 

G2 <- matrix(0, nrow = 6, ncol = 6)
G2[1:5,1:5] <- normG
G2[6,2] <- -normG[3,2]
G2[6,3] <- 1
D2 <- Lambda2 <- diag(6)
diag(D2)[1:5] <- diag(D)
diag(Lambda2)[1:5] <- diag(Lambda)

C2 <- solve(diag(6) - G2)%*%solve(D2)%*%Lambda2%*%solve(D2)%*%t(solve(diag(6) - G2))
C2[6,4]/C2[6,3]
normG[4,3]

max(abs(C2[1:5,1:5] - C))
-normG[3,2]*C[,2] + C[,3]
C2[,6]

# using Chen et al. 2014 to identify constraint (Theorem 1 and Lemma 1)

(C2[6,4] - G2[6,2]*C2[2,4] - G2[6,3]*C2[3,4])/(C2[6,3] - G2[6,2]*C2[2,3] - C2[6,3]*C2[3,3])
G2[4,3]


(diag(6) - G2)%*%C2

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

gg32 <- G[3,2]/(1-G[3,3]) # parameters known from dynamical observation

# id R23
ff <- function(x) 
  C[3,3] - (lambda[3]*((C[3,2] - gg32*(C[2,2] - lambda[3]*x^2))/(lambda[3]*x))^2 + gg32^2*(C[2,2] - lambda[3]*x^2))

plot(ff, 0, 10)

R23 <- uniroot(ff, c(0,1000), tol = .Machine$double.eps/1000)$root
ff(R[2,3])

R23
R[2,3]

# id R33

R33 <- (C[3,2] - gg32*(C[2,2] - lambda[3]*R23^2))/(lambda[3]*R23)
R33
R[3,3]

# id, what coef is this?
id1 <- (C[3,4] - gg32*C[4,2])/(lambda[3]*R33*(R33-gg32*R23))


###
# id R41/R31

# id R41

gg14 <- G[1,4]/(1-G[1,1]) # parameters known from dynamical observation

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

#
id1
id2

