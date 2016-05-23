# least squares adjusment mix model

## out function calculate the output for one iteration 
out <- function(B, A, L, P, W){
    ## funcion F(X, L) = 0. (modelo matematico)
    ## B, matriz (r,u), derivada de F con respecto a los parametros X
    ## A, matriz (r,n), derivada de F con respecto a las medidas L
    ## L, Mediciones vector (n,1)
    ## P, Matriz de pesos (n,n), (sigma^2)*inv(matriz VarCor)
    ## W, Vector (r,1) de la funcion F(X0, L0)
    require(MASS)
    QLb <- solve(P) # Matriz de Cofactores de las Observaciones tomadas
    r <- nrow(B) # Número de ecuaciones
    n <- ncol(B) # Número de variables
    M <- B %*% QLb %*% t(B) # Elemento ecuación normal
    N <- t(A) %*% solve(M) %*% A # Elemento ecuación normal
    U <- t(A) %*% solve(M) %*% W # Elemento ecuación normal
    X <- solve(-N) %*% U # Parámetros estimados
    NX <- -U # Ecuación normal
    K <- solve(-M) %*% (A %*% X + W) # Multiplicadores de Lagrange
    V <- QLb %*% t(B) %*% K # Residuos Estimados
    Qx <- solve(N) # Matriz de Cofactores de los parámetros
    Qv <- QLb %*% t(B) %*% solve(M) %*% 
          (M - A %*% solve(N) %*% t(A)) %*%
          solve(M) %*% B %*% QLb # Matriz de cofactores de los residuos
    QLa <- QLb - Qv # Matriz de Cofactores de las Observaciones Ajustadas
    
    
    return (list(X = c(X), V = c(V), L = L + c(V), Qv = Qv, QL = QLa, Qx = Qx))

}

## meval evaluates a matrix function output from sympy
## in the points peval, return a numeric matrix
meval <- function(m, peval){
         ## m es la matriz de expresiones salida de Rsympy
         ## peval es el vector de los valores apriori de arg
         ## arg = F(peval), vector caracteres de argumentos y mediciones apriori
         ## X(parámetros = x1, x2, ..., xn) y L0 (mediciones = l1, l2,..., ln)
         m <- as.matrix(m)    
         arg <- paste(names(peval), peval, sep = "=", collapse = ",")
         r <- nrow(m)
         c <- ncol(m)
         out <- matrix(data = rep(0.0, r*c), nrow = r, ncol = c)
         ## llenar la matrix de salida evaluando las funciones en los puntos peval
         for(r in seq_len(r)){
             
             for(c in seq_len(c)){
                 body <- m[r,c]
                 eval(parse(text = paste('f <- function(', arg, ') { return(' , body , ')}', sep='')))
                 out[r, c] <- f()
                 rm(f)
             }
             
         }
         return(out)
         
}

## mgen generates two matrix from a seed function and L, X vectors
## A derivates with respect to X parameters
## B derivates with respect to L measures

mgen <- function(chaFUN, L, X){
        # chaFUN a string vector with the generatrix Functions (F(X,L) = 0)
        # L vector de mediciones (l1, l2, ..., ln)
        # X vector de parámetros (x1, x2, ..., xn)
        require(rSymPy)
        L0 <- c(names(L))
        for(l0 in seq_along(L0)){
            eval(parse(text = paste(L0[l0], "<- Var('", L0[l0], "')", sep = "")))
        }
        X0 <- c(names(X))
        for(x0 in seq_along(X0)){
            eval(parse(text = paste(X0[x0], "<- Var('", X0[x0], "')", sep = "")))
        }
            
        nl <- length(L0)
        nx <- length(X0)
        r <- length(chaFUN)
        B <- matrix(as.numeric(rep(0.0, r*nl)), nrow = r, ncol = nl)
        A <- matrix(as.numeric(rep(0.0, r*nx)), nrow = r, ncol = nx)
        
        
        
        for(r in seq_len(r)){
            for(nl in seq_len(nl)){
                dl <- paste("diff(", chaFUN[r],",", L0[nl], ",1)")
                B[r,nl] <- sympy(dl)
            }
        }

        for(r in seq_len(r)){
            for(nx in seq_len(nx)){
                dx <- paste("diff(", chaFUN[r],",", X0[nx], ",1)")
                A[r,nx] <- sympy(dx)
            }
        }
        return(list(B = B, A = A))
        
}

## mix model adjusment
## generate an adjusment 

mmadjust <- function(chaFUN, L, X, VarCor, sigma2, maxiter = 10, tol = 1e-6){
    require(rSymPy)
    tol <- tol
    if(missing(sigma2) & missing(VarCor)){
        sigma2 <- 1
        VarCor <- diag(x = 1, nrow = length(L), ncol = length(L))
    }
    if(missing(sigma2) & !missing(VarCor)){
        sigma2 <- max(VarCor)
    }
    r <- length(chaFUN) # numero de filas
    P <- sigma2 * solve(VarCor)
    L0 <- L
    X0 <- X
    namesX <- names(X)
    namesL <- names(L)
    iter <- 1
    mgen0 <- mgen(chaFUN = chaFUN, L = L, X = X)
    B <- meval(mgen0$B, peval = c(L,X))
    A <- meval(mgen0$A, peval = c(L,X))
    W <- meval(m = chaFUN, peval = c(L,X))
    ajuste0 <- out(B = B, A = A, L = L, P = P, W = W)
    V0 <- ajuste0$L - L0
    r0 <- (t(V0)%*% P %*% V0)
    X <- ajuste0$X + X
    L <- ajuste0$L
    names(X) <- namesX
    names(L)<- namesL
    while(iter <= maxiter){
        mgen_i <- mgen(chaFUN = chaFUN, L = L, X = X)
        B <- meval(mgen_i$B, peval = c(L,X))
        A <- meval(mgen_i$A, peval = c(L,X))
        W <- meval(m = chaFUN, peval = c(L,X))
        ajuste_i <- out(B = B, A = A, L = L, P = P, W = W)
        V_i <- ajuste_i$L - L0
        r_i <- (t(V_i)%*% P %*% V_i)
        print("Residuo de Iteración: ", abs(r_i - r0))
        X <- ajuste_i$X + X
        L <- ajuste_i$L
        namesX <- names(X)
        namesL <- names(L)
        print(abs(r_i - r0))
        if(abs(r_i - r0) <= tol){
            message("Ajuste logrado en la iteración: ", iter)
            sigma2 <- r_i/(r - qr(A)$rank)
            break
        }
        if(abs(r_i - r0) > tol){
            r0 <- r_i
        }
        iter <- iter + 1
    }
    if(iter == maxiter + 1){
        message("No se logro Ajuste, iteración N°: ", iter)
        return(NULL)
    }
    if(iter < maxiter + 1){
        return(list(L = L, V = V_i, X = X,   
                    QL = ajuste_i$QL, 
                    Qx = ajuste_i$Qx,
                    Qv = ajuste_i$Qv,
                    sigma2 = sigma2,
                    W = W))
    }
}
    




