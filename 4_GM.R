dir="~"

require("ClusterR")
set.seed(1)


## ============
##  functions
## ============
gen.data <- function(n = 1000, cnt.rate = 0.01){
  n1 = ceiling(n * (1-cnt.rate))  
  n2 = n - n1
  
  z1 = rnorm(n=n1, mean=0, sd=1)
  z2 = rnorm(n=n1, mean=-5, sd=1)
  
  selection = rbinom(n=n1, size=1, prob=0.4)
  
  x1 = z1 * selection + z2 * (1-selection)
  x2 = rnorm(n=n2, mean=10, sd=0.1)
  
  return(list(n1=n1, n2=n2, x1=x1, x2=x2, 
              x=append(x1,x2), cnt.rate=cnt.rate))
}
pt <- function(theta, x){
  mu1 = theta[1]; sd1 = theta[2]
  mu2 = theta[3]; sd2 = theta[4]
  alpha = theta[5]
  alpha * dnorm(x=x, mean=mu1, sd=sd1) + (1-alpha) * dnorm(x=x, mean=mu2, sd=sd2)
}
tt <- function(theta,x){
  mu1 = theta[1]; sd1 = theta[2]
  mu2 = theta[3]; sd2 = theta[4]
  alpha = theta[5]
  coef1 = alpha * dnorm(x=x, mean=mu1, sd=sd1)/pt(theta,x)
  coef2 = (1-alpha) * dnorm(x=x, mean=mu2, sd=sd2)/pt(theta,x)
  v1 = coef1 * c((x-mu1)/(sd1^2), (x-mu1)^2/(sd1^3)-1/sd1)
  v2 = coef2 * c((x-mu2)/(sd2^2), (x-mu2)^2/(sd2^3)-1/sd2)
  c(v1,v2)
}
gen.y <- function(theta, m=10){
  mu1 = theta[1]; sd1 = theta[2]
  mu2 = theta[3]; sd2 = theta[4]
  alpha = theta[5]
  z1 = rnorm(n=m, mean=mu1, sd=sd1)
  z2 = rnorm(n=m, mean=mu2, sd=sd2)
  selection = rbinom(n=m, size=1, prob=alpha)
  z1 * selection + z2 * (1-selection)
}



mle <- function(x){
  gmm <- GMM(
    data=data.frame(x),
    gaussian_comps = 2,
    dist_mode = "eucl_dist",
    seed_mode = "random_subset",
    km_iter = 100,
    em_iter = 100,
    verbose = FALSE,
    var_floor = 1e-10,
    seed = 1,
    full_covariance_matrices = FALSE
  )
  
  mus = as.vector(gmm$centroids)
  sds = as.vector(sqrt(gmm$covariance_matrices))
  alpha = gmm$weights[1]
  c(mus[1], sds[1], mus[2], sds[2], alpha)
}




## =============
##  computation
## =============
cnt.rate = 0.01
data = gen.data(cnt.rate=cnt.rate)
save(file=paste0(dir,"/data_GM.RData"), data)

model = "GM"

beta_seq = c(0.5)
m_seq = c(10,25,100)
T=1000 #number of iterations
eta0=1; decay=0.7; decay_interval=25

settings = expand.grid(beta_seq, m_seq)

theta0 <- mle(x=data$x)

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







pp_mle = sapply(xx, function(x) pt(theta=theta0, x=x))
yl = range(0, list.pt, pp_mle)



pdf(file=paste0(dir,"/histogram_GM.pdf"), width=5, height=4)

par(mar = c(4, 4, 2, 4) + 0.3)  
hist(data$x, breaks=50, xlim=xl, main=" ", xlab="x")
par(new=T)
plot(xx, pp_mle, xlim=xl, type="l", xlab=" ",ylab=" ",xaxt="n", yaxt="n",
     col="black", lwd=2, ylim=yl)
colors =c(rep("blue",1), rep("red",1), rep("orange",1))
ltys = rep(c(1), 3)
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

