
# IUR calculation for SHR, SMR and STrR

## IUR for SHR

Data preparation:
```R
  obs_admission = "h_admits" #observed admission
  exp_admission = "expectda" #expeceted admission
  hosp_yar = "h_dy_yar" #hospital YAR
  facility = "provfs" #facility id
  year = "year" #data year
```
Code example:
```R

```


## Function lists

* Stratified sampling:
```R
  stra_sampling(fac)
    #sampling with replacement within each strata.
    #fac: a sorted factor vector, ex. facility.
```  
* IUR after bootsrap
```R
  IUR_bootdata(size,measure.boot,measure.org)
    #         size: facility size vector. Each item corresponds to a facility size.
    # measure.boot: bootstraped measure vector.
    #  measure.org: original measure vector, ex. SMR and SHR.
    
    # OUTPUT:
    #      IUR: total IUR
    #       nF: number of facilities
    # IUR.fac: faclity level IUR
```

* IUR_bootstrap
```R
  IUR_bootstrap(Obs, Exp, fac, n.boot = 100, stratify.var = NULL,stratify.cut = NULL, measure.fun = cal_SMR, seed = 123)
    # Note: IUR using bootstrap method. The data should be sorted by facitlity (fac).
    #          Obs: observed outcomes;
    #          Exp: expected outcomes;
    #          fac: facility id vector;
    #       n.boot: the number of bootstraps;
    # stratify.var: stratification variable;
    # stratify.cut: stratification cutoff. It should be a vector with two cutoff points. If it is NULL, tertiles will be used.
    #  measure.fun: function to calculate the measure of interest. For now, we only implement SMR type of function.
    #         seed: an integer number to generate random numbers.
```         
          
* cal_SMR
```R
  cal_SMR(Obs,Exp,fac)
    #Note: caculate SMR type measure. Data should be sorted by facitlity (fac).
    # Obs: observed outcomes;
    # Exp: expected outcomes;
    # fac: facility id vector;
```          
