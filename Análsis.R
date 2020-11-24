# Trabajo Final  ----------------------------------------------------------

##Distribución Gompertz Flexible Weibull (GoFW)
##(Se ajusta $\eta$ para linealizarla)
## $\eta = exp(\beta_0+\beta_1*x)$

library(dplyr)
library(readxl)
library(ggplot2)
library(stats4)
library(survival)


# Log - Verosimilitud -----------------------------------------------------

ll_gofw <- function(alpha, beta0, beta1, gamma, theta){
  n = length(y)
  vero = suppressWarnings(n*log(theta) + sum(log(alpha+(exp(beta0 + beta1*x))/y^2)) +
    sum(log(alpha*y^2-(exp(beta0 + beta1*x))/y)) +
    gamma * sum(exp(alpha*y-(exp(beta0 + beta1*x))/y)) + 
    (theta/gamma) * (1-(exp(-exp(alpha*y-(exp(beta0 + beta1*x))/y)))^-gamma))
  log_vero = -sum(vero)
  return(log_vero)
}

# Potencia ----------------------------------------------------------------


#Caso particular Parámetros igual a 1

gen_gofw <- function(n, b_0 = 2, b_1 = 0.5){
  u = runif(n)
  x =  rbinom(n,1,0.5)
  p = -log(log(1-log(1-u)))
  t = (-p+sqrt(p**2+4*exp(b_0+b_1*x)))/2
  return(list(t=t,x=x))
}


pot_gofw <- function(it=100, n){
  pot_b0 = NULL
  pot_b1 = NULL
  for (i in 1:it){
    datos = gen_gofw(n)
    m = mle(ll_gofw, start = list(alpha = 1, beta0 =0 , 
                                  beta1=0, gamma=1, theta= 1))
  pot_b0[i] = 1*(prod(summary(m)@coef[2,1]+c(-2,2)*summary(m)@coef[2,2]) > 0)
  pot_b1[i] = 1*(prod(summary(m)@coef[3,1]+c(-2,2)*summary(m)@coef[3,2]) > 0)
  }
  p_b0=mean(pot_b0)
  p_b1=mean(pot_b1)
  return(list(p_b0=p_b0, p_b1=p_b1))
}

n <- seq(1,500,10)
p1 <- NULL
p2 <- NULL
for(j in n){
  p1[j] <- pot_gofw(n=j)$p_b0
  p2[j] <- pot_gofw(n=j)$p_b1
}

# Selcción de modelo

data <- data.frame(
  t = gen_gofw(10)$t,
  e = c(rep(1,10)),
  x = gen_gofw(10)$x
)
so <- Surv(time = data$t, event = data$e)
m_wei <- survreg(so ~ data$x, data = so, dist="weibull")
m_exp <- survreg(so ~ data$x, data = so, dist="exponential")
y <- data$t
x <- data$x
m_gofw <- mle(ll_gofw, start = list(alpha = 1, beta0 =0 , 
                                    beta1=0, gamma=1, theta= 1))
c(AIC(m_gofw),AIC(m_wei),AIC(m_exp))

datos_pot <- read_excel("datos_pot.xlsx",sheet = "Hoja2")

datos_pot$beta <- factor(datos_pot$beta)
str(datos_pot)

ggplot(datos_pot,aes(x=n,y=potencia,group=beta)) +
  geom_line(aes(colour=beta),size=1.2,alpha=0.5) +
  geom_point(aes(colour=beta),size=2.2) +
  scale_color_brewer(palette="Dark2")+
  theme_bw()

avion <- read_excel("avion.xlsx", sheet = "Hoja1")

ggplot(avion, aes(x = t, color = as.factor(x),
                  fill = as.factor(x) )) +
  geom_histogram(aes(y =..density..), alpha=0.5, 
                 position="identity") +
  geom_density(alpha=.2) +
  xlab("Tiempo (años)") + ylab("Densidad") +
  scale_fill_discrete(name = "Tipo de Avión", 
                      labels = c("Tipo 1", "Tipo 2")) + 
  guides(color = FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal()

