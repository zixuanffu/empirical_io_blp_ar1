
## Standard
I always need to remind myself of the standard BLP estimation procedure.

The subscripts are $i$ for individual, $j$ for product, $t$ for time (market).

1. Given 
$$ u_{ijt}=x_{jt}\beta + \xi_{jt}+ \mu_{ijt} + \epsilon_{jt}=\delta_{jt}+\mu_{ijt}+\epsilon_{jt}$$
2. The random coefficient is contained in $\mu_{jt}$. Here I assume the only random coefficient is that of the constant term $\mu_{ijt}\sim N(0,\sigma^2)$. 

3. So we guess it first. Compute $\delta_{jt}$ from the observed market share $s_{jt}$.  
4. Regress $\delta_{jt}$ on $x_{jt}$ to get ${\beta}$ (with instruments).  
5. Having estimated the linear parameters $\beta$, we now move on to estimate the nonlinear one $\sigma$. We need to find another instrument for $\sigma$. We can construct GMM objective function written as $G(\sigma)WG(\sigma)'$.  
6. We minimize the GMM objective function to get an estimate of $\sigma$.

One can notice that the estimation of the linear parameters $\beta$ is done separately from the estimation of the nonlinear parameter $\sigma$. 
I would apply the same idea in the next case.

## AR(1) Structure

1. Given 
$$ u_{ijt}=x_{jt}\beta + \xi_{jt}+ \mu_{ijt} + \epsilon_{jt}=\delta_{jt}+\mu_{ijt}+\epsilon_{jt}$$
The unobserved product characteristics $\xi_{jt}$ follows an AR(1) process:
$$\xi_{jt} = \rho \xi_{j,t-1} + \omega_{jt}$$

2. Therefore, the non linear parameters are $\rho$ and $\sigma$. One from the AR(1) process and one from the random coefficient.

3. As before, we guess $\sigma$ first. Compute $\delta_{jt}$ from the observed market share $s_{jt}$.

4. Previously we have this moment condition
$$ E(x_{jt}(\delta_{jt}-x_{jt}\beta))=0$$ 
to estimate $\beta$. Now we have a new
equation since we have further imposed the AR(1) structure on the residual $\xi_{jt}$, which is
$$ \delta_{jt}-\rho \delta_{j,t-1} = \beta x_{jt} - \rho\beta x_{j,t-1} + \omega_{jt} $$
5. We first guess a value for $\rho$. We make use of this new equation which gives the moment condition to estimate the linear parameters $\beta$.
$$ E[x_{jt}(\delta_{jt}-\rho \delta_{j,t-1} - (\beta x_{jt} - \rho\beta x_{j,t-1}))]=0 $$
6. Having estimated the linear parameters $\beta$, we move on to estimating nonlinear ones, $\sigma$ and $\rho$. We need to find two additional instruments for them to construct the GMM objective function.
7. We minimize the GMM objective function $G(\sigma,\rho)WG(\sigma,\rho)'$ to get estimates of $\sigma$ and $\rho$.
8. Notice that the linear and nonlinear parameters are estimated separately. The last step is to estimate them jointly to get the correct standard errors.
We can construct $G(\beta,\sigma,\rho)WG(\beta,\sigma,\rho)'$ with the previous estimates (from previous separate estimation) as the initial values.


