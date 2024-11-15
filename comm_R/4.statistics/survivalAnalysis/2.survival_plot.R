install.packages(c("survival", "survminer"))


draw_survival_plot=function(indata,outcome,outcome_t,var){
  formula_in=formula(paste0("survival::Surv(",outcome_t,",",outcome,")","~",var))
  fit_one = do.call(survfit, args = list(formula = formula_in, data = indata))
  survminer::ggsurvplot(fit_one, data = indata, risk.table = F, pval = T)
}



#fit_one = survival::survfit(formula_in, data = indata)

#summary(fit_one)
#summary(fit_one)$table

#d <- data.frame(time = fit_one$time,
#                n.risk = fit_one$n.risk,
#                n.event = fit_one$n.event,
#                n.censor = fit_one$n.censor,
#                surv = fit_one$surv,
#                upper = fit_one$upper,
#                lower = fit_one$lower
#)




