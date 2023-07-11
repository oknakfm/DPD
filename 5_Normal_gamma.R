dir="~"

set.seed(2)

## ============
##  functions
## ============
gen.data <- function(n = 1000, cnt.rate = 0.01){
  n1 = ceiling(n * (1-cnt.rate))  
  n2 = n - n1
  
  x1 = rnorm(n=n1, mean=0, sd=1)
  x2 = rnorm(n=n2, mean=10, sd=0.1)
  
  return(list(n1=n1, n2=n2, x1=x1, x2=x2, 
              x=append(x1,x2), cnt.rate=cnt.rate))
}
pt <- function(theta, x){
  mu = theta[1]; sd = theta[2]
  dnorm(x=x, mean=mu, sd=sd)
}
tt <- function(theta,x){
  mu = theta[1]; sd = theta[2]
  c((x-mu)/(sd^2), (x-mu)^2/(sd^3)-1/sd)
}
gen.y <- function(theta, m=10){
  mu = theta[1]; sd = theta[2]
  rnorm(n=m, mean=mu, sd=sd)
}


## =============
##  computation
## =============
cnt.rate = 0.1
data = gen.data(cnt.rate=cnt.rate)

save(file=paste0(dir,"/data_normal.RData"), data)

gamma_seq = c(0.1,0.3,0.5)
m_seq = c(3,10,100)
T=500 #number of iterations
eta0=1; decay=0.7; decay_interval=25

theta0 <- c(mean(data$x), sd(data$x), 1)
settings = expand.grid("normal", gamma_seq, m_seq)

for(setting_id in 1:nrow(settings)){
  model = settings[setting_id,1]
  gamma = settings[setting_id,2]
  m = settings[setting_id,3]
  
  ## optimization
  theta = theta0
  recorded.theta <- theta0

  mu = theta[1]; sd = theta[2]
  pt_X = sapply(data$x, function(x) pt(theta, x))
  rt = (2*pi*sd^2)^(-gamma/2) * (1+gamma)^(-1/2)
  loss = -log(mean(pt_X^gamma))/gamma + log(rt)/(1+gamma)
  
  recorded.loss = loss
  
  eta = eta0
  
  for(itr in 1:T){
    Y = gen.y(theta, m=m)
    
    pt_X = sapply(data$x, function(x) pt(theta, x))
    tt_X = t(sapply(data$x, function(x) tt(theta, x)))
    
    pt_Y = sapply(Y, function(y) pt(theta, y))
    tt_Y = t(sapply(Y, function(y) tt(theta, y)))
    
    scale = theta[3]
    stoc_grad = append(-scale^gamma * apply(pt_X^gamma * tt_X, 2, mean) + scale^(1+gamma) * apply(pt_Y^gamma * tt_Y, 2, mean),
                       -scale^(gamma-1) * mean(pt_X^gamma) + scale^gamma * mean(pt_Y^gamma))
    
    if(itr%%decay_interval==0) eta = eta * decay
    theta <- theta - eta * stoc_grad
    
    recorded.theta <- rbind(recorded.theta, theta)

    mu = theta[1]; sd = theta[2]
    rt = (2*pi*sd^2)^(-gamma/2) * (1+gamma)^(-1/2)
    loss = -log(mean(pt_X^gamma))/gamma + log(rt)/(1+gamma)
    recorded.loss = append(recorded.loss, loss)
  }
  name = paste0("normal_g_",gamma,"_",m)
  save(file=paste0(dir,"/",name,".RData"), recorded.theta, recorded.loss, theta)  
}


## ==================
##  plot 1 (fitting)
## ==================

xl = range(data$x)
xx = seq(xl[1], xl[2], length.out=200)

list.theta = NULL; list.pt = NULL
settings = expand.grid(gamma_seq, m_seq)
for(setting_id in 1:nrow(settings)){
  gamma = settings[setting_id,1]
  m = settings[setting_id,2]
  name = paste0("normal_g_",gamma,"_",m)
  load(file=paste0(dir,"/",name,".RData"))
  
  list.theta <- rbind(list.theta, theta)
  pp = sapply(xx, function(x) pt(theta=theta, x=x))
  list.pt <- rbind(list.pt, pp)
}
pp_mle = sapply(xx, function(x) pt(theta=c(mean(data$x), sd(data$x), 1), x=x))
yl = range(0, list.pt, pp_mle)

pdf(file=paste0(dir,"/histogram_normal_gamma.pdf"), width=9, height=8)

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

names = "MLE (gamma=0)"
for(setting_id in 1:nrow(settings)){
  names = append(names, paste0("gamma=",settings[setting_id,1],
                               " / m=",settings[setting_id,2]))
}

legend("topright", legend=names, 
       col = append("black", colors),
       lty = append(1, ltys))

dev.off()


## ===============
##  plot 2 (loss)
## ===============

for(gamma in gamma_seq){
  pdf(file=paste0(dir,"/GCE_",gamma,".pdf"), width=4, height=4)
  par(mar = c(4, 4, 2, 1), oma = c(0,0,0,0)) 
  list.loss = NULL
  for(m in m_seq){
    name = paste0("normal_g_",gamma,"_",m)
    load(file=paste0(dir,"/",name,".RData"))
    list.loss = rbind(list.loss, recorded.loss)
  }
  itr = 0:T
  xl = range(itr)
  yl = range(list.loss)
  
  colors=c("blue","red","orange")
  
  
  plot(0,mean(yl),xlim=xl, ylim=yl, xlab="t", ylab="GCE",
       main=paste0("gamma=",gamma), type="n")
  for(m_id in 1:length(m_seq)){
    par(new=T)
    plot(itr, list.loss[m_id,], type="l", xlim=xl, ylim=yl,
         xaxt="n",yaxt="n",xlab=" ",ylab=" ",
         col=colors[m_id])
  }
  legend("topright", legend=paste0("m=",m_seq), col=colors, 
         lty=1)
  
  dev.off()
}


## ===============================================
for(gamma in gamma_seq){
  pdf(file=paste0(dir,"/Gamma-scale_",gamma,".pdf"), width=4, height=4)
  par(mar = c(4, 4, 2, 1), oma = c(0,0,0,0)) 
  list.scale = NULL
  for(m in m_seq){
    name = paste0("normal_g_",gamma,"_",m)
    load(file=paste0(dir,"/",name,".RData"))
    list.scale = rbind(list.scale, recorded.theta[,3])
  }
  itr = 0:T
  xl = range(itr)
  yl = range(list.scale, 1-cnt.rate)
  
  colors=c("blue","red","orange")
  
  plot(0,mean(yl),xlim=xl, ylim=yl, xlab="t", ylab="scale parameter c>0",
       main=paste0("gamma=",gamma), type="n")
  for(m_id in 1:length(m_seq)){
    par(new=T)
    plot(itr, list.scale[m_id,], type="l", xlim=xl, ylim=yl,
         xaxt="n",yaxt="n",xlab=" ",ylab=" ",
         col=colors[m_id])
  }
  abline(h=1-cnt.rate)
  legend("topright", legend=paste0("m=",m_seq), col=colors, 
         lty=1)
  
  dev.off()
}





