L <- c(t1 = 0.6435, t2 = 0.9827, t3 = 1.5152, d1 = 13**0.5*(1.01), d2 = 5*(1.02) , d3 = 6*(1.03))
X <- c(x1 = 0.01, y1 = 0.025, x2 = 6.05, y2 = -0.1, x3 = 3.89, y3 = 2.93)
VarCor <- diag(c(rep(2.424068e-03, 3), (0.01 + 0.01*L[4:6])), ncol = 6)
sigma2 <- max(VarCor)
chaFUN <- c(f1 = "(t1 + t2 +t3) + (atan((y2 - y1)/(x2 - x1)) - atan((y3 - y1)/(x3 - x1)) + atan((y3 - y2)/(x3 - x2))- atan((y1 - y2)/(x1 - x2) + atan((y1 - y3)/(x1 - x3)) - atan((y2 - y3)/(x2 - x3)))) - 2*pi",
            f2 = "d3**2-(d1**2+d2**2-2*d1*d2*cos(t1))") 

mmadjust(chaFUN = chaFUN, L = L, X = X, VarCor = VarCor, sigma2 = sigma2)
gen <- mgen(chaFUN, L, X)
evalgenB <- meval(gen$B, peval = c(L,X))
evalgenA <- meval(gen$A, peval = c(L,X))