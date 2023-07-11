  dir="~"
  
  require("VGAM")
  set.seed(1)
  
  
  ## ============
  ##  functions
  ## ============
  gen.data <- function(n = 1000, cnt.rate = 0.01){
    n1 = ceiling(n * (1-cnt.rate))  
    n2 = n - n1
    
    x1 = rgompertz(n=n1, scale=1, shape=0.1)
    x2 = rnorm(n=n2, mean=10, sd=0.1)
    
    return(list(n1=n1, n2=n2, x1=x1, x2=x2, 
                x=append(x1,x2), cnt.rate=cnt.rate))
  }
pt <- function(theta, x){
  scale = theta[1]; shape = theta[2]
  dgompertz(x=x, scale=scale, shape=shape)
}
tt <- function(theta,x){
  scale = theta[1]; shape = theta[2]
  coef = (1-exp(scale*x))/(scale^2) + x*exp(scale*x)/scale
  c(x-shape*coef, 1/shape + (1-exp(scale*x))/scale)
}
gen.y <- function(theta, m=10){
  scale = theta[1]; shape = theta[2]
  rgompertz(n=m, scale=scale, shape=shape)
}


mle <- function(x){
  f <- function(w){
    weights = 1/mean(1-exp(w*x))
    v1 = (1-exp(w*x))/w
    v2 = x * exp(w*x)
    sum(x) + weights * (sum(v1) + sum(v2))
  }
  scale = uniroot(f, c(0.01,10))$root
  shape = -scale/mean(1-exp(scale*x))
  return(c(scale, shape))
}



## =============
##  computation
## =============
cnt.rate = 0.01
data = gen.data(cnt.rate=cnt.rate)

save(file=paste0(dir,"/data_gompertz.RData"), data)

model = "gompertz"

beta_seq = c(1)
m_seq = c(10, 25, 100)
T=1000 #number of iterations
eta0=0.5; decay=0.7; decay_interval=25

settings = expand.grid(beta_seq, m_seq)
theta0 = mle(x=data$x)


for(setting_id in 1:nrow(settings)){
  beta = settings[setting_id,1]
  m = settings[setting_id,2]
  
  ## optimization
  mu = mean(data$x)
  theta = theta0
  recorded.theta <- theta0
  
  pt_X = sapply(data$x, function(x) pt(theta, x))
  eta = eta0
  
  for(itr in 1:T){
    Y = gen.y(theta, m=m)
    
    pt_X = sapply(data$x, function(x) pt(theta, x))
    tt_X = t(sapply(data$x, function(x) tt(theta, x)))
    
    pt_Y = sapply(Y, function(y) pt(theta, y))
    tt_Y = t(sapply(Y, function(y) tt(theta, y)))
    
    stoc_grad = - apply(pt_X^beta * tt_X, 2, mean) + apply(pt_Y^beta * tt_Y, 2, mean)
    
    if(itr%%decay_interval==0) eta = eta * decay
    theta <- abs(theta - eta * stoc_grad)
    
    recorded.theta <- rbind(recorded.theta, theta)
    
  }
  name = paste0(model,"_",beta,"_",m)
  save(file=paste0(dir,"/",name,".RData"), recorded.theta, theta)  
}


## ==================
##  plot 1 (fitting)
## ==================

xl = range(data$x)
xx = seq(xl[1], xl[2], length.out=200)

list.theta = NULL; list.pt = NULL
settings = expand.grid(beta_seq, m_seq)
for(setting_id in 1:nrow(settings)){
  beta = settings[setting_id,1]
  m = settings[setting_id,2]
  name = paste0(model,"_",beta,"_",m)
  load(file=paste0(dir,"/",name,".RData"))
  
  list.theta <- rbind(list.theta, theta)
  pp = sapply(xx, function(x) pt(theta=theta, x=x))
  list.pt <- rbind(list.pt, pp)
}

pp_mle = sapply(xx, function(x) pt(theta=theta0, x=x))
yl = range(0, list.pt, pp_mle)

pdf(file=paste0(dir,"/histogram_gompertz.pdf"), width=5, height=4)

par(mar = c(4, 4, 2, 4) + 0.3)  
hist(data$x, breaks=50, xlim=xl, main=" ", xlab="x")
par(new=T)
plot(xx, pp_mle, xlim=xl, type="l", xlab=" ",ylab=" ",xaxt="n", yaxt="n",
     col="black", lwd=2, ylim=yl)
colors =c(rep("blue",1), rep("red",1), rep("orange",1))
ltys = rep(1, 3)
for(setting_id in 1:nrow(settings)){
  par(new=T)
  plot(xx, list.pt[setting_id,], xlim=xl, type="l", col=colors[setting_id],
       xlab=" ", ylab=" ", xaxt="n", yaxt="n", lty=ltys[setting_id], lwd=2, 
       ylim=yl)
}

axis(side = 4, at = pretty(yl))      # Add second axis
mtext("Probability density", side = 4, line = 3)  

names = "MLE (beta=0)"
for(setting_id in 1:nrow(settings)){
  names = append(names, paste0("beta=",settings[setting_id,1],
                               " / m=",settings[setting_id,2]))
}

legend("topright", legend=names, 
       col = append("black", colors),
       lty = append(1, ltys))

dev.off()

