b <- rcauchy(1000, location = 0, scale = 2.5)
c <- b[b>=0]
N <- length(c)

table(round(c,0))
plot(density(c))

a=rnorm(length(c),0,10)
plot(density(a))
table(round(a,0))

plot(density(a*c))

dc <- NULL
for (i in 1:length(N)){
  dc[i]=rnorm(1,0,1/c[i])

}

plot(density(dc))
mean(a*c)
mean(dc)


g <- rgamma(N,0.001,0.001)
dg <- NULL
for (i in 1:length(N)){
  dg[i]=rnorm(1,0,1/g)

}