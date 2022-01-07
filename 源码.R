## 数据与函数定义：

n<-c(55,157,159,16)
x<-c(0,2,7,3)
z<-matrix(c(1,1,1,1,7,14,27,57),nrow=4)
# 似然函数
logP=function(theta){
  result=0
  for(i in 1:4){
    result=result+x[i]*z[i,]%*%t(theta)-n[i]*log(1+exp(z[i,]%*%t(theta)))
    return(result)
  }
}
# Jacobi向量
J=function(theta){
  result=c(0,0)
  for(i in 1:4){
    for(j in 1:2){
      result[j]=result[j]+(x[i]-n[i]*exp(z[i,]%*%theta)/(1+exp(z[i,]%*%theta)))*z[i,j]
    }
  }
  return(result)
}
# Hession矩阵
H=function(theta){
  result=matrix(rep(0,4),nrow=2)
  for(i in 1:4){
    for(j in 1:2){
      for(k in 1:2){
        tmp=z[i,j]*z[i,k]
        tmpnum=n[i]*exp(z[i,]%*%theta)/(1+exp(z[i,]%*%theta))^2
        result[j,k]=result[j,k]-tmpnum[1,1]*tmp
      }
    }
    
  }
  return(result)
}


## 1. Newton's method:


# Newton法求根
findRoot=function(theta){
  delta=1000
  count=0
  while(sum(delta^2)>1e-6){
    Hmatrix=H(theta)
    Jmatrix=J(theta)
    delta=solve(Hmatrix)%*%Jmatrix
    theta=theta-0.1*delta  # step size = 0.1
    count=count+1
  }
  # print(count)  # 迭代次数
  return(theta)
}
findRoot(c(-5,0))

## 2. BFGS法：
findRoot2=function(theta){
  delta=1000
  count=0
  B=matrix(c(-11,-338,-338,-13322),nrow=2)
  while(sum(delta^2)>1e-6){
    delta=0.1*solve(B)%*%J(theta)
    thetap=theta-0.1*solve(B)%*%J(theta)
    s=thetap-theta
    y=J(thetap)-J(theta)
    if((t(y)%*%s)[1,1]!=0){
      B=B-B%*%s%*%t(s)%*%B/((t(s)%*%B%*%s)[1,1])+y%*%t(y)/((t(y)%*%s)[1,1])
    }
    theta=thetap  # 更新theta
    count=count+1
  }
  print(count)
  return(theta)
}
findRoot2(c(-5,0))


## 3. GLM package:

library(glmnet)

glmr=data.frame(
  z=c(7,14,27,57),n=c(55,157,159,16),x=c(0,2,7,3)
)
glmr$y<-cbind(glmr$x,glmr$n-glmr$x)
glm.solve<-glm(y~z,family = binomial,data=glmr)
summary(glm.solve)
