# Project: IncVAUCAP

This project provides the R code for the numerical study in the preprint <https://arxiv.org/abs/2010.09822>. 

## R packages required

`pracma`, `nleqslv`, `dplyr`, `parallel`, `doSNOW`, `foreach`

## Numerical Study

In this study, we evaluate the increment value (IncV) of adding a marker, denoted by <img src="https://render.githubusercontent.com/render/math?math=Y">, to a model with an existing marker, denoted by <img src="https://render.githubusercontent.com/render/math?math=X">. The markers <img src="https://render.githubusercontent.com/render/math?math=X"> and <img src="https://render.githubusercontent.com/render/math?math=Y"> be independent standard normal random variables. Given the values of these two markers, a binary outcome <img src="https://render.githubusercontent.com/render/math?math=D"> follows a Bernoulli distribution with the probability of <img src="https://render.githubusercontent.com/render/math?math=D=1"> via the following model: 

![](https://latex.codecogs.com/gif.latex?%5Cpi%28X%2CY%29%20%3D%20Pr%28D%3D1%5Cmid%20X%2CY%29%20%3D%20%5CPhi%28%5Cbeta_0&plus;%5Cbeta_1X&plus;%5Cbeta_2Y&plus;%5Cbeta_3XY%29)

where<img src="https://render.githubusercontent.com/render/math?math=\Phi(\cdot)"> is the CDF of a standard normal distribution. Given <img src="https://render.githubusercontent.com/render/math?math=X"> and <img src="https://render.githubusercontent.com/render/math?math=Y">, <img src="https://render.githubusercontent.com/render/math?math=\pi(X,Y)"> is the *true* risk. 

Typically, in practice, none of the working models are the true model. Having this in mind, we compare the following two misspecified working models: (i) \textit{one-marker model}: <img src="https://render.githubusercontent.com/render/math?math=p(X) = \Phi(\gamma_0+\gamma_1X)">, and (ii) \textit{two-marker model}: <img src="https://render.githubusercontent.com/render/math?math=p(X,Y) = \Phi(\gamma_0+\gamma_1X+\gamma_2Y)">. 

In this project, we consider three accuracy measures: AUC (area under the receiver operating characteristic curve), AP (area under the precision-recall curve), and sBrS (scaled Brier score). Their IncV parameters are $\Delta \Psi = \Psi_{M_2} - \Psi_{M_1}$, where $\Psi$ is AUC, AP, or sBrS, and $\Psi_{M_2}$ and $\Psi_{M_1}$ are the accuracy measures for the two-marker model and one-marker model respectively.  

Note that we are interested in the IncV parameters of the population working risk, not in the IncV estimates from a sample, we do not use simulation studies, and thus, the IncV parameters can be directly derived from the distributional assumptions described above.

## Files

* `helpers.R` includes 

    + `R` function which calculates the $\beta_0$ given the values of $\beta_1$, $\beta_2$, $\beta_3$, and the event rate $\pi_1=Pr(D=1)$;
    + `R` functions which calculates the values of the parameters $\gamma_0$, $\gamma_1$, and $\gamma_2$ for the one-marker and two-marker models;
    + `R` functions which calucates the distributional functions (such as cumulative distribution function, probability density function, and survival function) of the one-marker working risk score $r(X)=\gamma_0 + \gamma_1 X$ and two-marker working risk score $r(X,Y)=\gamma_0 + \gamma_1 X+\gamma_2Y$ for (i) the whole population, (ii) events (subjects with $D=1$), and (iii) non-events (subjects with $D=0$);
    + `R` functions which calculates AUC, AP, and Brier score.
    
* `numstudy.R` includes an `R` function which calculates the AUC, AP, and sBrS of the one-marker and two-marker models as well as their IncV, given the values of $\pi_1$, $\beta_1$, $\beta_2$, and $\beta_3$. For example 

```{r}
source("numstudy.R")
numstudy(pi1=0.01, beta1=1, beta2=0.8, beta3=0.2)
```
* `example.R` includes the R code for calculating the AUC, AP, and sBrS of the prescribed-dose model and ovarian-dose model (for predicting acute ovarian failure) as well as the IncV of ovarian-dose model, compared to the prescribed-dose model. In addition, the data is also provided in `dataAOF.csv`.
