#' test_wilcox
#'
#' @param df_in 
#' @param var 
#' @param g 
#' @param label 
#'
#' @return
#' @export
#'
#' @examples
test_wilcox=function(df_in,var,g,label="wilcox"){
  df_in1=df_in[c(var,g)]
  names(df_in1)=c("var","g")
  test_w=wilcox.test(var~g,data=df_in1)
  data.frame(label=label,var=var,g=g,statistic_wilcox=test_w$statistic,P_wilcox=test_w$p.value,P_wilcox_sci=signif(test_w$p.value,2))
}












