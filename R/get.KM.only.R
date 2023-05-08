
## Kaplan-Meier (without influnce function)
get.KM.only<-function(Xi,Di, t0.all){
  Ni.t  = outer(Xi,t0.all,FUN = "<=")*Di
  wei   = Ni.t+outer(Xi,t0.all,FUN = ">")
  G_fit = survfit(Surv(Xi,1-Di)~1,type="kaplan-meier")
  G_hat = t(apply(outer(Xi,t0.all,FUN=pmin),1,
                  function(x) summary(G_fit,times=x,extend = TRUE)$surv))
  KM    = apply(wei/G_hat*outer(Xi,t0.all,FUN=">"),2,sum)/apply(wei/G_hat,2,sum)

  return(KM)
}
