# Project: IncVAUCAP

This project provides the R code for the numerical study in the preprint <https://arxiv.org/abs/2010.09822>. 

## R packages required

`pracma`, `nleqslv`, `dplyr`, `parallel`, `doSNOW`, `foreach`

## Numerical Study

In this study, we evaluate the increment value (IncV) of adding a marker, denoted by <img src="https://render.githubusercontent.com/render/math?math=Y">, to a model with an existing marker, denoted by <img src="https://render.githubusercontent.com/render/math?math=X">. The markers <img src="https://render.githubusercontent.com/render/math?math=X"> and <img src="https://render.githubusercontent.com/render/math?math=Y"> be independent standard normal random variables. Given the values of these two markers, a binary outcome <img src="https://render.githubusercontent.com/render/math?math=D"> follows a Bernoulli distribution with the probability of <img src="https://render.githubusercontent.com/render/math?math=D=1"> via the following model: 

![](https://latex.codecogs.com/gif.latex?%5Cpi%28X%2CY%29%20%3D%20Pr%28D%3D1%5Cmid%20X%2CY%29%20%3D%20%5CPhi%28%5Cbeta_0&plus;%5Cbeta_1X&plus;%5Cbeta_2Y&plus;%5Cbeta_3XY%29)

where <img src="https://render.githubusercontent.com/render/math?math=\Phi(\cdot)"> is the CDF of a standard normal distribution. Given <img src="https://render.githubusercontent.com/render/math?math=X"> and <img src="https://render.githubusercontent.com/render/math?math=Y">, <img src="https://render.githubusercontent.com/render/math?math=\pi(X,Y)"> is the *true* risk. 

Typically, in practice, none of the working models are the true model. Having this in mind, we compare the following two misspecified working models: (i) *one-marker model*: ![](https://latex.codecogs.com/gif.latex?p%28X%29%20%3D%20%5CPhi%28%5Cgamma_0%20&plus;%20%5Cgamma_1X%29), and (ii) *two-marker model*: ![](https://latex.codecogs.com/gif.latex?p%28X%2CY%29%20%3D%20%5CPhi%28%5Cgamma_0&plus;%5Cgamma_1X&plus;%5Cgamma_2Y%29). 

In this project, we consider three accuracy measures: AUC (area under the receiver operating characteristic curve), AP (area under the precision-recall curve), and sBrS (scaled Brier score). Their IncV parameters are <img src="https://render.githubusercontent.com/render/math?math=\Delta \Psi = \Psi_{M_2} - \Psi_{M_1}">, where <img src="https://render.githubusercontent.com/render/math?math=\Psi"> is AUC, AP, or sBrS, and <img src="https://render.githubusercontent.com/render/math?math=\Psi_{M_2}"> and <img src="https://render.githubusercontent.com/render/math?math=\Psi_{M_1}"> are the accuracy measures for the two-marker model and one-marker model respectively.  

Note that we are interested in the IncV parameters of the population working risk, not in the IncV estimates from a sample, we do not use simulation studies, and thus, the IncV parameters can be directly derived from the distributional assumptions described above.

## Files

* `helpers.R` includes 

    + `R` function which calculates the <img src="https://render.githubusercontent.com/render/math?math=\beta_0"> given the values of <img src="https://render.githubusercontent.com/render/math?math=\beta_1, \beta_2, \beta_3">, and the event rate <img src="https://render.githubusercontent.com/render/math?math=\pi_1 = Pr(D=1)">;
    + `R` functions which calculates the values of the parameters <img src="https://render.githubusercontent.com/render/math?math=\gamma_0, \gamma_1, \gamma_2"> for the one-marker and two-marker models;
    + `R` functions which calucates the distributional functions (such as cumulative distribution function, probability density function, and survival function) of the one-marker working risk score ![](https://latex.codecogs.com/gif.latex?r%28X%29%20%3D%20%5Cgamma_0&plus;%5Cgamma_1X) and two-marker working risk score ![](https://latex.codecogs.com/gif.latex?r%28X%2CY%29%20%3D%20%5Cgamma_0&plus;%5Cgamma_1X&plus;%5Cgamma_2Y) for (i) the whole population, (ii) events (subjects with <img src="https://render.githubusercontent.com/render/math?math=D=1">), and (iii) non-events (subjects with <img src="https://render.githubusercontent.com/render/math?math=D=0">);
    + `R` functions which calculates AUC, AP, and Brier score.
    
* `numstudy.R` includes an `R` function which calculates the AUC, AP, and sBrS of the one-marker and two-marker models as well as their IncV, given the values of <img src="https://render.githubusercontent.com/render/math?math=\pi_1, \beta_1, \beta_2, \beta_3">. For example 

```{r}
source("numstudy.R")
numstudy(pi1=0.01, beta1=1, beta2=0.8, beta3=0.2)
```
* `AOF_example.R` includes the R code for estimating the AUC, AP, and sBrS of the prescribed-dose model and ovarian-dose model (for predicting acute ovarian failure (AOF)) as well as generating the plots in Figure 1 of the manuscript. In addition, the data is also provided in `AOF.csv`.
