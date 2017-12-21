# Math 324 cheatsheet, Fall 2017 McGill, By Andre Kaba

* **Linear estimation**:

    $S_{xx}=\left[\sum_{i=1}^{n}x_i^2\right]-n(\overline{x}_n)^2,S_{yy}=\left[\sum_{i=1}^{n}y_i^2\right]-n(\overline{y}_n)^2,S_{xy}=\left[\sum_{i=1}^{n}x_iy_i\right]-n(\overline{x}_n\overline{y}_n)$

    $\hat{\beta}_1=\frac{S_{xy}}{S_{xx}},\hat{\beta_0}=\overline{y}_n-\hat{\beta_1}\overline{x}_n,\hat{\sigma}^2=\frac{S_{yy}-(\hat{\beta}_1)^2S_{xx}}{n-2},Var(\beta_0)=\hat{Var(\hat{\beta_1})}\left(\frac{\sum_{i=1}^{n} x_i^2}{n}\right), \hat{Var(\hat{\beta}_1)}=\frac{\hat{\sigma}^2}{S_{xx}}$

* **Hyp Testing**:
    Common tests: (dont forget to say under $H_0$ when writing $\sim$)

    * Single:
        * $x_1,...x_n \sim N(\mu,\sigma) \text{ where the other param is unknown }$
            * $\mu \text{ vs } \mu_0$, $n\geq 30 \implies \frac{\overline{X}_n-\mu_0}{S_n/\sqrt{n}}\approx N(0,1) \text{ otherwise } \frac{\overline{X}_n-\mu_0}{S_n/\sqrt{n}}\sim t_{(n-1)}$
            * $\sigma \text{ vs } \sigma_0$, $\frac{(n-1)S^2_n}{\sigma^2_0} \sim \chi^2_{(n-1)}$
        * $x_1,...x_n \sim Bern(p)$
            * $p \text{ vs } p_0$, $n \geq 30, \implies \frac{\hat{p}-p_0}{\sqrt{\frac{p_0(1-p_0)}{n}}}\approx N(0,1)$

    * Multi (dont forget to mention indep $X_i$,$Y_i$)
         * $X_1,...X_n \sim N(\mu_0,\sigma),Y_1,...Y_m \sim N(\mu_1,\sigma) \text{ where the other param is unknown }$
            * $\mu_0-\mu_1 \text{ vs } 0$, $n\geq 30 \implies \frac{\overline{X}_n-\overline{Y}_m}{\sqrt{\frac{S^2_n}{n}+\frac{S^2_m}{m}}}\approx N(0,1) \text{ otherwise } \frac{\overline{X}_n-\overline{Y}_m}{S_{pooled}\sqrt{\frac{1}{n}+\frac{1}{m}}}\sim t_{(n+m-2)},S^2_{pooled}=\frac{(m-1)S_m^2+(n-1)S^2_n}{m+n-2}$
        * $X_1,...X_n \sim Bern(p_0),Y_1,...Y_m \sim Bern(p_1)$
            * $p_0-p_1 \text{ vs } 0$, $n \geq 30, \implies \frac{\hat{p}_0-\hat{p}_1}{\sqrt{\hat{p}(1-\hat{p})(\frac{1}{m}+\frac{1}{n})}}\approx N(0,1),\hat{p}=\frac{X+Y}{n+m}$

    Np test for only one missing variable: $RR:\{\frac{L(\theta_1)}{L(\theta_0)}>k\}$

    Likelihood ration (LR) test for one or more missing variables: $RR:\{\frac{max_{\theta \in \Theta_0}L(\theta)}{max_{\theta \in \Theta}L(\theta)}=\lambda(X) \leq C_{\alpha}\}, C_{\alpha} \in [0,1]$

    With regularity conditions for large $n$:

    $2\left(max_{\theta \in \Theta}[ln(L(\theta)]-max_{\theta \in \Theta_0}[ln(L(\theta)]\right)=-2ln(\lambda(X))$

    $\implies RR:\{ -2ln[\lambda(X)] \geq C^*_{\alpha}=-2ln(C_{\alpha})\}\approx \{ \chi^2_d \sim \chi^2_{(d)} \geq C^*_{\alpha}\}, d=dim(\theta)-\dim(\theta_0)$

    **If we cannot reject we say**: based on the given data, we do not have enough evidence to reject the null hypothesis at the significance level $\alpha$

* **Conf Intervals**:
    * If the interval is concrete i.e $[a,b]$ with $95\%$ we say:
    * If we repeat this experiment 100 times we expect $95\%$ of these intervals to contain $\theta$
    * The prexperiment CI that lead to these values has $95\%$ chance of containing $\theta$
    
    * $Uni(0,\theta) \implies \left[\frac{X_{(n)}}{\theta}\right]^n \implies \theta \in \left[\frac{X_{(n)}}{(1-\alpha/2)^{1/n}},\frac{X_{(n)}}{(\alpha/2)^{1/n}}\right]$

    * $N(\mu,\sigma) \implies \mu \in \left[\overline{X}_n\pm t_{(n-1);\alpha/2}\frac{S_n}{\sqrt{n}}\right], \sigma^2 \in \left[\frac{(n-1)S_n^2}{\chi^2_{(n-1);\alpha/2}},\frac{(n-1)S_n^2}{\chi^2_{(n-1);1-\alpha/2}}\right], \mu_0-\mu_1 \in \left[ (\overline{X}_n-\overline{Y}_m) \pm t_{(n+m-2);\alpha/2}S_{pooled}\sqrt{\frac{1}{n}+\frac{1}{m}}\right]$

    * Large samples: $\mu \in \mu \in \left[\overline{X}_n\pm z_{\alpha/2}\frac{S_n}{\sqrt{n}}\right], \mu_0-\mu_1 \in \left[ (\overline{X}_n-\overline{Y}_m) \pm z_{\alpha/2}\sqrt{\frac{S_n^2}{n}+\frac{S_m^2}{m}}\right], p \in \left[\hat{p}_n \pm z_{\alpha/2}\sqrt{\frac{\hat{p}_n(1-\hat{p}_n)}{n}}\right]$
    $p_0-p_1 \in \left[ (\hat{p}_0-\hat{p}_1) \pm z_{\alpha/2}\sqrt{\frac{\hat{p}_0(1-\hat{p}_0)}{n}+\frac{\hat{p}_1(1-\hat{p}_1)}{m}}\right], \theta \in \left[\hat{\theta}\pm z_{\alpha/2}\sqrt{\frac{1}{n}[\hat{I(\theta)}]^-1}\right]$

* **MLE Methods**
    * Moments:
        * $E(X^k)=m_k=\frac{1}{n}\sum_{i=1}^{n} X_i^k$, we need $k$ equations for $k$ parameters
    * MLE:
        * we compute $\hat{\theta}_n=argmax_{\theta \in \Theta} ln[L(\theta)]$, it's often a function of MSS by factorization theorm and if it's biased we fix it to obtain a MVUE
        * Invariance prop: if $\hat{\theta}_n$ is MLE for $\theta$ then $\tau(\hat{\theta}_n)$ is a MLE of $\tau(\theta)$
        * regualar conditions and large sample: $\sqrt{n}(\hat{\theta}_n-\theta) \xrightarrow{d} N(0,I^{-1}(\theta))$, $I^{-1}(\theta)=-E\left[\frac{\partial^2 \ln f(X;\theta)}{\partial \theta^2}\right]=-\sum_{i=1}^{n} \frac{\partial^2}{\partial \theta^2}\ln f(x_i;\theta)$ is fisher info matrix

* **Estimation (Don't forget the $\mathbb{1}_{Y_{(n)}\leq \theta_1}\mathbb{1}_{Y_{(1)}\geq \theta_0}$ for Uni)**:
    * $MSE(\hat{\theta}_n)=E[(\hat{\theta}_n-\theta)^2]=Var(\hat{\theta}_n)+[B(\hat{\theta}_n)]^2,eff(\hat{\theta}_n,\hat{\theta'}_n)=\frac{Var(\hat{\theta}'_n)}{Var(\hat{\theta}_n)}$
    * we want consistant estimators. i.e $\lim_{n \rightarrow \infty} P(|\hat{\theta}_n-\theta)\leq \epsilon)=1$ for any $\epsilon$ or $\lim_{n \rightarrow \infty} Var(\hat{\theta})=0$. This means that $\hat{\theta}_n \xrightarrow{p} \theta$ as $n \rightarrow \infty$
    * $T$ is sufficient statistic for $\theta$ if $f(X_1,...,X_n|T)$ doesn't depend on $\theta$, or by **Factorization Theorem** $t=T(X_1,)$ is sufficient $\iff L(\theta)=g(t,\theta)h(x_1,...,x_n),g$ only depend on $x_1,...,x_n$ through $T$, $h$ doesn't depend on $\theta$
    * **Minimal sufficient statistics** using Lehmann Criterion: $T$ is minimal sufficient if:
        For any two points $(x_1,...,x_n),(y_1,...,y_n), \frac{L(x_1,...,x_n)}{L(y_1,...,y_n)}$ doesn't depend on $\theta$ iff $T(x_1,...,x_n)=T(y_1,...,y_n)$
        * Normal $\mu \implies \sum_{i=1}^{n} x_i,\sigma^2 \implies (\sum_{i=1}^{n} x_i,\sum_{i=1}^{n} x_i^2)$/Uni $(\theta_1,\theta_1) \implies (Y_{(1)},Y_{(n)})$
        * Poisson/Binomial $\lambda \implies \sum_{i=1}^{n} x_i, p \implies \sum_{i=1}^{n} x_i$
        * Gamma $\alpha \implies \sum_{i=1}^{n} x_i, \beta \implies \prod_{i=1}^{n} x_i$/ Beta $\alpha \implies \prod_{i=1}^{n} x_i, \beta \implies \prod_{i=1}^{n} (1-x_i)$

    * Rao-blackwell, we can use a MSS $T$ to improve an unbiased estimator $\hat{\theta}$ as follows, $\hat{\theta'}=E[\hat{\theta}|T]$, it will have same bias but less or equal variance, tzkr $E(f(Y)|Y=y)=\sum_{i=1}^{n} \frac{P(f(Y)=i,Y=y)}{P(Y=y)}$
    * $\hat{\theta}$ is MVUE if it's unbiased and a function of a sufficient statistic (mo darori minimal) $T$, $\hat{\theta}=g(T)$
    * Normal approx of binomial: $P(X\geq x)\approx P(Z \geq x-0.5),P(X\leq x)\approx P(Z \leq x+0.5)$ (continuity correction), always do $P(Y>y)=1-P(Y\leq y)\approx 1-P(X\leq y+0.5), X \sim N(np,np(1-p))$
    * Estimation with bounds $P(|\hat{\theta}-\theta|\leq 2\sigma^2_{\hat{\theta}})=0.95$ and we say "the probability that the error of this estimation within the 2-standard error is $2\sigma^2_{\hat{\theta}}$"

* **Distros**:
    * $X \sim Exp(\beta) \implies P(X\geq b)=e^{-b/\beta}$
    * $\chi^2_{r_1}+...+\chi^2_{r_n}=\chi^2_{r_1+...+r_n}$
    * $Gamma(\alpha_1,\beta)+...+Gamma(\alpha_n,\beta)=Gamma(\alpha_1+...+\alpha_n,\beta), Exp(\beta)=Gamma(1,\beta),\chi^2_{(r)}=Gamma(r/2,2)$
    * $Poi(\lambda_1)+...+Poi(\lambda_n)=Poi(\lambda_1+...+\lambda_n)$)
    * $\Gamma(\alpha)=\int_{0}^{\infty} x^{\alpha-1}e^{-x} dx$,$\Gamma(\alpha+1)=\alpha \Gamma(\alpha), \forall \alpha>0$, $\Gamma(n)=(n-1)!$,$\Gamma(1/2)=\sqrt{\pi}$
    * $B(\alpha,\beta)=\int_0^{1} x^{\alpha-1}(1-x)^{\beta-1} dx$,$B(\alpha,\beta)=\frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha+\beta)}$
    * Distro-----mean-----variance
    * Binomial: $P(k)={n \choose k}p^{k}(1-p)^{n-k},np,np(1-p)$/Poisson: $P(y)=\frac{\lambda^{y}e^{-\lambda}}{y!},\lambda,\lambda$/Uniform: $f(y)=\frac{1}{\theta_2-\theta_1}\mathbb{1}_{\theta_1 \leq y \leq \theta_2},\frac{\theta_1+\theta_2}{2},\frac{(\theta_2-\theta_1)^2}{12}$
    * Normal: $f(y)=\frac{1}{\sigma\sqrt{2\pi}}exp\left[-\left(\frac{1}{2\sigma^2}\right)(y-\mu)^2 \right]$/Gamma: $f(y)=\left[\frac{1}{\Gamma(\alpha)\beta^{\alpha}}\right]y^{\alpha-1}e^{-y/\beta},\alpha\beta,\alpha\beta^2$/Exponential: $f(y)=\frac{1}{\beta}e^{-y/\beta};\beta>0,\beta,\beta^2$/Geometric $p(y)=p(1-p)^{y-1},\frac{1}{p},\frac{1-p}{p^2}$