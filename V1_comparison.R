require("mvtnorm","RColorBrewer")
dir = "~"
d = 2

n = 500
beta = 0.5
dev = 100
cnt.rate = 0.01
theta.true = rep(0.5,d)

## for SGD
T_SGD = 300 #number of iterations
eta0=0.5; decay=0.7; decay_interval=20
m_seq = c(3, 10, 50)

## for Fullbatch GD
T_FGD = 300
omega0 = mean(eta0 * decay^(0:((T_SGD/decay_interval)-1)))
M_seq = c(3, 10, 50)
D = 2 ## NI is performed over [-D,D]^d



## ============
## functions
## ============
pt <- function(theta, X){
  n = nrow(X)
  p = sapply(1:n, function(i)
    dmvnorm(x=X[i,], mean=theta, sigma=diag(d))
  )
  return(p)
}
tt <- function(theta,X){
  t(t(X)-theta)
}
gen.y <- function(theta, m=10){
  rmvnorm(n=m, mean=theta, sigma=diag(d))
}

n_pur = ceiling(n * (1-cnt.rate))
n_cnt = n - n_pur

X_pur = rmvnorm(n=n_pur, mean=theta.true, sigma=diag(d))
X_cnt = rmvnorm(n=n_cnt, mean=theta.true + rep(dev,d), sigma=0.01*diag(d))
X = rbind(X_pur,X_cnt)


theta0 = apply(X,2,mean)

## SGD
for(m in m_seq){
  theta = theta0
  eta = eta0
  error = NULL
  for(itr in 1:T_SGD){
    Y = gen.y(theta=theta, m=m)
    stoc_grad = -apply(tt(theta,X) * pt(theta,X)^beta, 2, mean) + apply(tt(theta,Y) * pt(theta,Y)^beta, 2, mean)
    
    if(itr%%decay_interval==0) eta = eta * decay
    theta <- theta - eta * stoc_grad
    error = append(error, sum((theta - theta.true)^2))
  }
  complexity = (n+m) * (1:T_SGD)
  name = paste0(dir,"/",d,"_SGD_",m,".RData")
  save(file=name, complexity, error)
}

## fullbatch GD
for(M in M_seq){
  xx = seq(-D,D,length.out=M)
  Y = expand.grid(data.frame(matrix(rep(xx,d), M, d)))
  theta = theta0
  omega = omega0
  error = NULL
  for(itr in 1:T_FGD){
    Y = gen.y(theta=theta, m=m)
    grad = -apply(tt(theta,X) * pt(theta,X)^beta, 2, mean) + apply((2*D)^d * tt(theta,Y) * pt(theta,Y)^beta, 2, mean)
    theta <- theta - omega * grad
    error = append(error, sum((theta - theta.true)^2))
  }
  complexity = (n+(M^d)) * (1:T_FGD)
  name = paste0(dir,"/",d,"_FGD_",M,".RData")
  save(file=name, complexity, error)
}



SGD_complexities = SGD_errors = NULL
for(m in m_seq){
  name = paste0(dir,"/",d,"_SGD_",m,".RData")
  load(file=name)
  SGD_complexities = rbind(SGD_complexities, complexity)
  SGD_errors = rbind(SGD_errors, error)
}

FGD_complexities = FGD_errors = NULL
for(M in M_seq){
  name = paste0(dir,"/",d,"_FGD_",M,".RData")
  load(file=name)
  FGD_complexities = rbind(FGD_complexities, complexity)
  FGD_errors = rbind(FGD_errors, error)
}

xl = range(SGD_complexities,FGD_complexities)
xl[2] = xl[2] * 10
yl = range(SGD_errors,FGD_errors)
logtype="x"


pdf(file=paste0(dir,"/d=",d,".pdf"), width=5, height=4)
par(cex = 1)
par(mar = c(4, 4, 1, 1), oma = c(0,0,0,0)) 

plot(mean(xl),mean(yl), xlim=xl, ylim=yl,
     xlab="complexity", ylab="error", type="n",
     log=logtype)

n_SGD = nrow(SGD_complexities)
SGD_colors = brewer.pal(n=3+n_SGD, name="BuPu")[3+(1:n_SGD)]

for(i in 1:n_SGD){
  par(new=T)
  plot(SGD_complexities[i,],SGD_errors[i,], xlim=xl, ylim=yl,
       xlab=" ", ylab=" ", xaxt="n", yaxt="n", type="l", log=logtype,
       col = SGD_colors[i])
}

n_FGD = nrow(FGD_complexities)
FGD_colors = brewer.pal(n=3+n_FGD, name="YlGn")[3+(1:n_FGD)]
for(i in 1:n_FGD){
  par(new=T)
  plot(FGD_complexities[i,],FGD_errors[i,], xlim=xl, ylim=yl,
       xlab=" ", ylab=" ", xaxt="n", yaxt="n", type="l", log=logtype,
       col = FGD_colors[i], lty=5)
}
par(cex = 0.75)
legend("topright", legend=c(paste0("SGD / m=",m_seq),
                            paste0("GD / M=",M_seq^d)),
       lty=c(rep(1,n_SGD),rep(5,n_FGD)),
       col=c(SGD_colors,FGD_colors), bg="white")
dev.off()