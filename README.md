# logistic-regression-
统计计算课程大作业。用R语言实现。

# 一、问题描述：

拟合 quantal response 模型。quantal response 模型的covariate vector为$z$，系数为$\theta$，观测变量$y_i$服从二项分布。具体来说，第i个观测变量为成功的概率为：
$$
P(y_i=1;\theta)=\pi_i(\theta)=\frac{e^{z_i^T\theta}}{1+e^{z_i^T\theta}}
$$
此即为logistic模型。可以发现这是广义线性模型：

连接函数为：

$$
z_i^T\theta=g(\mu_i)=log\frac{\mu_i}{1-\mu_i}
$$
且$y_i$属于指数分布族：

$$
P[y_i=1|z_i]=\frac{e^{z_i^T\theta}}{1+e^{z_i^T\theta}}
$$

# 二、数据：

对同一个$z_i$做$n_i$次实验，不妨记观测变量为$\{y_{ij}\},j=1..n_i$。令$x_i=\sum_{j}y_{ij}$。于是我们的数据表示如下：

实验次数n=(55,157,159,16),成功次数x=(0,2,7,3),z=(7,14,27,57)。
  
# 三、公式计算：

根据分布函数，容易计算出以下结果。

## 1. 必要函数：
### a. 概率密度函数

$$
P(z,x,n;\theta)=\prod_{i=1}^4 \frac{\left(e^{z_i^T\theta}\right)^{x_i}}{\left({1+e^{z_i^T\theta}}\right)^{n_i}}
$$

### b. 似然函数

$$
l(z,x,n;\theta)=logP=\sum_{i=1}^4 \left(x_iz_i^T\theta -n_ilog(1+e^{z_i^T\theta})\right)
$$

### c. Jacobi向量

$$
J=\sum_{i=1}^4 \left(x_i -\frac{n_ie^{z_i^T\theta}}{1+e^{z_i^T\theta}}\right)z_i
$$

### d. Hessian 矩阵

$$
H=
-\sum_{i=1}^4\frac{n_ie^{z_i^T\theta}}{(1+e^{z_i^T\theta})^2} z_{i}\times z_{i}^T
$$

## 2. 优化方法：

为了求解最优的$\theta$，我们想找到令$l$最大化的$\hat{\theta}$。下面我们分别使用三种优化方法求解此优化问题。

### a. Newton法：

$$
\theta^{(t+1)}
=
\theta^{(t)}
-H(\theta^{(t)})^{-1}
J(\theta^{(t)})
$$

### b. Fisher Scoring法：

$$
\theta^{(t+1)}
=
\theta^{(t)}
+I(\theta^{(t)})^{-1}
J(\theta^{(t)})
$$

其中 I的计算如下，

$$
I(\theta)=-\mathbb{E}[H|\theta]=-H
$$
可以发现：针对本问题，由于$H$不含$x$，求期望前后形式一致。因此此问题Newton法与Scoring法的迭代公式一致。

### c. Quasi-Newton法：

Quasi-Newton法的思想是近似$B\approx H$，有多种实现方式。这里我们不妨采用常见的方法之一，BFGS。算法如下：

$$
\theta^{(t+1)}=\theta^{(t)}-B(\theta^{(t)})^{-1}J(\theta^{(t)})
$$
然后更新$B$：

$$
B(\theta^{(t+1)})=B(\theta^{(t)})-\frac{B(\theta^{(t)})s_ts_t^TB(\theta^{(t)})}{s_t^TB(\theta^{(t)})s_t}+\frac{y_ty_t^T}{y_t^Ts_t}
$$

其中
$$
s_t=\theta^{(t+1)}-\theta^{(t)},y_t=J(\theta^{(t+1)})-J(\theta^{(t)})
$$


# 四、编程：

采用上述算法，其中迭代停止的标志为$\Delta\theta$的$l2$范数小于$\epsilon$。具体见源码。
