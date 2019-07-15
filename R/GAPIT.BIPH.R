#Object: To generate binary phenotype
 #Straitegy: NA
 #Output: binary phenotype (0 and 1's)
 #intput: genetic effect (x), hertiability (h2) and ratio of 1's (r)
 #Authors: Zhiwu Zhang
 #Last update: March 18, 2016
##############################################################################################
`GAPIT.BIPH` <-
function(x=0,h2=.5,r=.25){
    #To assign probability for given standard normal variable x and h2
    #Author: Zhiwu Zhang
    #Last update: Febuary 27, 2016
    p=pnorm(x)
    srp=1-p-r
    sh=1/(1-sqrt(h2))
    adj=(r-.5)*(1-sqrt(h2))
    f=1/(1+exp(sh*srp))+adj
    return(f)
  }
#=============================================================================================


