mylink <- function(eta,link){
  switch(link,
         "probit" = {return(pnorm(eta))},
         "logit" = {return(exp(eta)/(1+exp(eta)))}
  )
}
beta0_find = function(pi1,beta,link="probit"){
  fun.obj <- function(beta0,link){
    fun0 = function(x,y,beta0,link){
      eta0 = beta0 + beta[1]*x + beta[2]*y + beta[3]*x*y
      pi0 = mylink(eta0,link=link)
      pi0*dnorm(x)*dnorm(y)
    }
    pracma::integral2(fun0,xmin=-8,xmax=8, ymin=-8, ymax=8, beta0=beta0,link=link)$Q - pi1
  }
  uniroot(fun.obj,lower=-10,upper=10,link=link)$root
}


ParaEst_OneMarker <- function(beta_true,link="probit"){
  
  glm.fun0 = function(x1,x2,beta_true,gamma,index,link){
    eta_true = beta_true[1] + beta_true[2]*x1 + beta_true[3]*x2 + beta_true[4]*x1*x2
    eta_fit = gamma[1] + gamma[2]*x1; 
    pi_true = mylink(eta_true,link=link)
    pi_fit = mylink(eta_fit,link=link)
    switch(link,
           "probit" = {out=dnorm(eta_fit)*(1/(pi_fit*(1-pi_fit)))*(pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)},
           "logit" = {out = (pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)}
    )
    return(out)
  }
  
  est.fun <- function(gamma,beta_true,link=link){
    s0 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(0,0), link=link)$Q
    s1 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(1,0),link=link)$Q
    c(s0,s1)
  }
  nleqslv::nleqslv(c(0,0), est.fun, beta_true=beta_true,link=link)$x
}

ParaEst_TwoMarker <- function(beta_true,link="probit"){
  
  glm.fun0 = function(x1,x2,beta_true,gamma,index,link){
    eta_true = beta_true[1] + beta_true[2]*x1 + beta_true[3]*x2 + beta_true[4]*x1*x2
    eta_fit = gamma[1] + gamma[2]*x1 + gamma[3]*x2; 
    pi_true = mylink(eta_true,link=link)
    pi_fit = mylink(eta_fit,link=link)
    switch(link,
           "probit" = {out=dnorm(eta_fit)*(1/(pi_fit*(1-pi_fit)))*(pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)},
           "logit" = {out = (pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)}
    )
    return(out)
  }
  
  est.fun <- function(gamma,beta_true, xmin,xmax,ymin,ymax,link){
    s0 = pracma::integral2(glm.fun0,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                           gamma=gamma,beta_true = beta_true, index=c(0,0), link=link)$Q
    s1 = pracma::integral2(glm.fun0,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                           gamma=gamma,beta_true = beta_true, index=c(1,0), link=link)$Q
    s2 = pracma::integral2(glm.fun0,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                           gamma=gamma,beta_true = beta_true, index=c(0,1), link=link)$Q
    c(s0,s1,s2)
  }
  cc = 8; index = T
  while (index & cc>=4){
    out = try(nleqslv::nleqslv(c(0,0,0), est.fun, beta_true=beta_true, xmin=-cc,xmax=cc,ymin=-cc,ymax=cc,link=link),silent = TRUE)
    index = is.character(out)
    cc = cc - 1
  }
  if (index) return(rep(NA,3)) else return(out$x)
}

ParaEst_TrueModel <- function(beta_true,link="probit"){
  
  glm.fun0 = function(x1,x2,beta_true,gamma,index,link){
    eta_true = beta_true[1] + beta_true[2]*x1 + beta_true[3]*x2 + beta_true[4]*x1*x2
    eta_fit = gamma[1] + gamma[2]*x1 + gamma[3]*x2 + gamma[4]*x1*x2; 
    pi_true = mylink(eta_true,link=link)
    pi_fit = mylink(eta_fit,link=link)
    switch(link,
           "probit" = {out=dnorm(eta_fit)*(1/(pi_fit*(1-pi_fit)))*(pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)},
           "logit" = {out = (pi_true - pi_fit)*x1^index[1]*x2^index[2]*dnorm(x1)*dnorm(x2)}
    )
    return(out)
  }
  
  est.fun <- function(gamma,beta_true){
    s0 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(0,0))$Q
    s1 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(1,0))$Q
    s2 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(0,1))$Q
    s3 = pracma::integral2(glm.fun0,xmin=-7,xmax=7,ymin=-7,ymax=7,
                           gamma=gamma,beta_true = beta_true, index=c(1,1))$Q
    c(s0,s1,s2,s3)
  }
  nleqslv::nleqslv(c(0,0,0,0), est.fun, beta_true=beta_true)$x
}

yy_pdf <- function(yy,gamma,given="x2"){
  ## Pr(Y=y)
  sapply(yy,function(y){
    if (given=="x2"){
      Jointpdf <- function(x2,y,gamma){
        ## Pr(X2=x2 & Y=y)
        x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
        (1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
      }
      pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
    } else {
      Jointpdf <- function(x1,y,gamma){
        ## Pr(X2=x2 & Y=y)
        x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
        (1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
      }
      pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
    }
  })
}


yy_pdf_case <- function(yy, gamma, beta_true, link="probit", given="x2"){
  ## Pr(D=1 & Y=y)
  sapply(yy,function(y){
    if (given=="x2"){
      Jointpdf_case <- function(x2,y,gamma,beta_true,link){
        ## Pr(D=1 & X2=x2 & Y=y)
        x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
        eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
        out = mylink(eta,link=link)*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
        out[is.na(out)] = 0
        return(out)
      }
      pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
    } else{
      Jointpdf_case <- function(x1,y,gamma,beta_true,link){
        ## Pr(D=1 & X2=x2 & Y=y)
        x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
        eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
        out = mylink(eta,link=link)*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
        out[is.na(out)] = 0
        return(out)
      }
      pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
    }
  })
}

yy_pdf_control <- function(yy, gamma, beta_true, link="probit", given="x2"){
  ## Pr(D=0 & P=p)
  sapply(yy,function(y){
    if (given=="x2"){
      Jointpdf_control <- function(x2,y,gamma,beta_true,link){
        ## Pr(D=0 & X2=x2 & Y=y)
        x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
        eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
        out = (1-mylink(eta,link=link))*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
        out[is.na(out)]=0
        return(out)
      }
      pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
    } else {
      Jointpdf_control <- function(x1,y,gamma,beta_true,link){
        ## Pr(D=0 & X2=x2 & Y=y)
        x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
        eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
        out = (1-mylink(eta,link=link))*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
        out[is.na(out)]=0
        return(out)
      }
      pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
    }
    
  })
}

yy_cdf = function(yy,gamma,y.min,given="x2"){
  ## Pr(Y<=y)
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    yy_pdf <- function(yy,gamma,given="x2"){
      ## Pr(Y=y)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf <- function(x2,y,gamma){
            ## Pr(X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            (1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
          }
          pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
        } else {
          Jointpdf <- function(x1,y,gamma){
            ## Pr(X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            (1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
          }
          pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
        }
      })
    }
    pracma::integral(yy_pdf,y.min,yy[k],gamma=gamma, given=given)
  }
  return(out)
}

yy_survival = function(yy,gamma,y.max, given="x2"){
  ## Pr(Y>=p)
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    yy_pdf <- function(yy,gamma,given="x2"){
      ## Pr(Y=y)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf <- function(x2,y,gamma){
            ## Pr(X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            (1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
          }
          pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
        } else {
          Jointpdf <- function(x1,y,gamma){
            ## Pr(X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            (1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
          }
          pracma::integral(Jointpdf, -8, 8,y=y,gamma=gamma)
        }
      })
    }
    pracma::integral(yy_pdf,yy[k],y.max,gamma=gamma, given=given)
  }
  return(out)
}


yy_cdf_case = function(yy,gamma,beta_true, y.min, link="probit", given="x2"){
  ## Pr(D=1 & Y <= y )
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    mylink <- function(eta,link){
      switch(link,
             "probit" = {return(pnorm(eta))},
             "logit" = {return(exp(eta)/(1+exp(eta)))}
      )
    }
    yy_pdf_case <- function(yy, gamma, beta_true, link="probit", given="x2"){
      ## Pr(D=1 & Y=y)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf_case <- function(x2,y,gamma,beta_true,link){
            ## Pr(D=1 & X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = mylink(eta,link=link)*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
            out[is.na(out)] = 0
            return(out)
          }
          pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
        } else{
          Jointpdf_case <- function(x1,y,gamma,beta_true,link){
            ## Pr(D=1 & X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = mylink(eta,link=link)*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
            out[is.na(out)] = 0
            return(out)
          }
          pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
        }
      })
    }
    
    pracma::integral(yy_pdf_case,y.min,yy[k],gamma=gamma,beta_true=beta_true,link=link, given=given)
  }
  return(out)
}
yy_survival_case = function(yy,gamma,beta_true,y.max, link="probit", given="x2"){
  ## Pr(D=1 & Y >= y )
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    mylink <- function(eta,link){
      switch(link,
             "probit" = {return(pnorm(eta))},
             "logit" = {return(exp(eta)/(1+exp(eta)))}
      )
    }
    yy_pdf_case <- function(yy, gamma, beta_true, link="probit", given="x2"){
      ## Pr(D=1 & Y=y)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf_case <- function(x2,y,gamma,beta_true,link){
            ## Pr(D=1 & X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = mylink(eta,link=link)*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
            out[is.na(out)] = 0
            return(out)
          }
          pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
        } else{
          Jointpdf_case <- function(x1,y,gamma,beta_true,link){
            ## Pr(D=1 & X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = mylink(eta,link=link)*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
            out[is.na(out)] = 0
            return(out)
          }
          pracma::integral(Jointpdf_case,-5, 5,y=y,gamma=gamma, beta_true=beta_true, link=link)
        }
      })
    }
    
    pracma::integral(yy_pdf_case,yy[k],y.max,gamma=gamma,beta_true=beta_true,link=link, given=given)
  }
  return(out)
}

yy_cdf_control = function(yy,gamma,beta_true, y.min, link="probit", given="x2"){
  ## Pr(D=0 & Y <= y )
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    mylink <- function(eta,link){
      switch(link,
             "probit" = {return(pnorm(eta))},
             "logit" = {return(exp(eta)/(1+exp(eta)))}
      )
    }
    yy_pdf_control <- function(yy, gamma, beta_true, link="probit", given="x2"){
      ## Pr(D=0 & P=p)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf_control <- function(x2,y,gamma,beta_true,link){
            ## Pr(D=0 & X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = (1-mylink(eta,link=link))*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
            out[is.na(out)]=0
            return(out)
          }
          pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
        } else {
          Jointpdf_control <- function(x1,y,gamma,beta_true,link){
            ## Pr(D=0 & X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = (1-mylink(eta,link=link))*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
            out[is.na(out)]=0
            return(out)
          }
          pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
        }
        
      })
    }
    
    pracma::integral(yy_pdf_control,y.min,yy[k],gamma=gamma,beta_true=beta_true, link=link, given=given)
  }
  return(out)
}
yy_survival_control = function(yy,gamma,beta_true, y.max, link="probit", given="x2"){
  ## Pr(D=0 & Y <= y )
  out = foreach(k=1:length(yy),.combine=c,.packages = "pracma") %dopar% {
    mylink <- function(eta,link){
      switch(link,
             "probit" = {return(pnorm(eta))},
             "logit" = {return(exp(eta)/(1+exp(eta)))}
      )
    }
    yy_pdf_control <- function(yy, gamma, beta_true, link="probit", given="x2"){
      ## Pr(D=0 & P=p)
      sapply(yy,function(y){
        if (given=="x2"){
          Jointpdf_control <- function(x2,y,gamma,beta_true,link){
            ## Pr(D=0 & X2=x2 & Y=y)
            x1 = (y-gamma[3]*x2-gamma[1])/(gamma[2]+gamma[4]*x2)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = (1-mylink(eta,link=link))*(1/abs(gamma[2]+gamma[4]*x2))*dnorm(x1)*dnorm(x2)
            out[is.na(out)]=0
            return(out)
          }
          pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
        } else {
          Jointpdf_control <- function(x1,y,gamma,beta_true,link){
            ## Pr(D=0 & X2=x2 & Y=y)
            x2 = (y-gamma[2]*x1-gamma[1])/(gamma[3]+gamma[4]*x1)
            eta = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
            out = (1-mylink(eta,link=link))*(1/abs(gamma[3]+gamma[4]*x1))*dnorm(x1)*dnorm(x2)
            out[is.na(out)]=0
            return(out)
          }
          pracma::integral(Jointpdf_control,-8, 8,y=y,gamma=gamma, beta_true = beta_true,link=link)
        }
        
      })
    }
    
    pracma::integral(yy_pdf_control,yy[k],y.max,gamma=gamma,beta_true=beta_true, link=link, given=given)
  }
  return(out)
}


auc_fun = function(gamma,pi1,beta_true,y.min,y.max, link="probit", given="x2"){
  pi0 = 1-pi1
  auc_obj_fun = function(yy,gamma,link, given){
    return(yy_cdf_control(yy,gamma,beta_true,y.min, link=link, given=given)*yy_pdf_case(yy,gamma,beta_true, link=link, given=given))
  }
  return(pracma::integral(auc_obj_fun,y.min,y.max,gamma=gamma,link=link, given=given)/(pi1*pi0))
}

ap_fun = function(gamma,pi1,beta_true,y.min,y.max, link="probit", given="x2"){
  ap_obj_fun = function(yy,gamma,link, given){
    return((yy_survival_case(yy,gamma,beta_true,y.max, link=link, given=given)/yy_survival(yy,gamma,y.max, given=given))*(yy_pdf_case(yy,gamma,beta_true, link=link, given=given)/pi1))
  }
  return(pracma::integral(ap_obj_fun,y.min,y.max,gamma=gamma, link=link, given=given))
}

bs_fun = function(gamma,beta_true,link="probit"){
  fun0 = function(x1,x2,link){
    eta_true = beta_true[1]+beta_true[2]*x1+beta_true[3]*x2+beta_true[4]*x1*x2
    pi_true = mylink(eta_true,link=link)
    eta_fit = gamma[1] + gamma[2]*x1 + gamma[3]*x2
    pi_fit = mylink(eta_fit,link=link)
    return((pi_true*(1-pi_true)+(pi_true-pi_fit)^2)*dnorm(x1)*dnorm(x2))
  }
  bs = pracma::integral2(fun0,-8,8,-8,8,link=link)$Q
  return(bs)
}

