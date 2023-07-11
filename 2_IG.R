dir="~"

require("actuar")
set.seed(1)


## ============
##  functions
## ============
gen.data <- function(n = 1000, cnt.rate = 0.01){
  n1 = ceiling(n * (1-cnt.rate))  
  n2 = n - n1
  
  x1 = rinvgauss(n=n1, mean=1, shape=3)
  x2 = rnorm(n=n2, mean=10, sd=0.1)
  
  return(list(n1=n1, n2=n2, x1=x1, x2=x2, 
              x=append(x1,x2), cnt.rate=cnt.rate))
}
pt <- function(theta, x){
  mu = theta[1]; shape = theta[2]
  dinvgauss(x=x, mean=mu, shape=shape)
}
tt <- function(theta,x){
  mu = theta[1]; shape = theta[2]
  c(shape*(x-mu)/(mu^3), 1/(2*shape) - (x-mu)^2/(2*mu^2*x))
}
gen.y <- function(theta, m=10){
  mu = theta[1]; shape = theta[2]
  rinvgauss(n=m, mean=mu, shape=shape)
}


## =============
##  computation
## =============
cnt.rate = 0.1
data = gen.data(cnt.rate=cnt.rate)

save(file=paste0(dir,"/data_IG.RData"), data)

model = "IG"

beta_seq = c(0.1,0.5,1)
m_seq = c(3,10,100)
T=1000 #number of iterations
eta0=1; decay=0.7; decay_interval=25

settings = expand.grid(beta_seq, m_seq)

theta0 <- c(mu, 1/mean(1/data$x - 1/mu)) ## MLE

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
    theta <- theta - eta * stoc_grad
    
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

mu = mean(data$x)
theta0 <- c(mu, 1/mean(1/data$x - 1/mu))
pp_mle = sapply(xx, function(x) pt(theta=theta0, x=x))
yl = range(0, list.pt, pp_mle)

pdf(file=paste0(dir,"/histogram_IG.pdf"), width=9, height=8)

par(mar = c(4, 4, 2, 4) + 0.3)  
hist(data$x, breaks=50, xlim=xl, main=" ", xlab="x")
par(new=T)
plot(xx, pp_mle, xlim=xl, type="l", xlab=" ",ylab=" ",xaxt="n", yaxt="n",
     col="black", lwd=2, ylim=yl)
colors =c(rep("blue",3), rep("red",3), rep("orange",3))
ltys = rep(c(1,2,4), 3)
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

