######################################## Cochran-Armitage test for trend, One group variable(Oridinal, more than 2 levels), One catgory variable ########################################
cal_CochranArmitageTest_p=function(indata,var_group,var_x_list){
  bind_rows(lapply(var_x_list, function(var_cat_one){
    x=table(indata[c(var_cat_one,var_group)])
    
    if(nrow(x)==1){
      out=data.frame(
        stringsAsFactors = F,
        Variable=var_cat_one,
        variable_levels=row.names(x)[1],
        test_name=NA,
        P_trend_raw=NA
      )
    }
    
    
    if(nrow(x)==2){
      test_out=DescTools::CochranArmitageTest(x, "two.sided")
      out=data.frame(
        stringsAsFactors = F,
        Variable=rep(var_cat_one,2),
        variable_levels=row.names(x)[1:2],
        test_name=rep(test_out$method,2),
        P_trend_raw=rep(signif(test_out$p.value,2),2)
      )
      
    }
    
    if(nrow(x)>=3){
      out=bind_rows(
        lapply(1:nrow(x), function(i){
          x1=x[i,]
          x2=apply(x[-i,],2,sum)
          x3=rbind(x1,x2)
          test_out=DescTools::CochranArmitageTest(x3, "two.sided")
          out=data.frame(
            stringsAsFactors = F,
            Variable=var_cat_one,
            variable_levels=row.names(x)[i],
            test_name=test_out$method,
            P_trend_raw=signif(test_out$p.value,2)
          )
          out
        })
      )
    }
    out$P_trend=convert_p(out$P_trend_raw)
    out
  }))
  
}
