library(forestmodel)
library(survival)
library(dplyr)
library(visR)
pretty_lung <- lung %>%
  transmute(time,
            status,
            Age = age,
            Sex = factor(sex, labels = c("Male", "Female")),
            ECOG = factor(lung$ph.ecog)) 

coxph(Surv(time, status) ~ ., pretty_lung)
print(forest_model(coxph(Surv(time, status) ~ ., pretty_lung)))


lung_cohort <- survival::lung %>%
  dplyr::transmute(AVAL = time,
                   Status = status - 1,
                   Sex = factor(sex, labels = c("Male", "Female")),
                   ECOG = factor(survival::lung$ph.ecog)) 

visR_HR <- lung_cohort %>%
  visR::estimate_KM(strata = c("Sex", "ECOG"), CNSR = "Status") %>%
  my_func(update_formula = "survival::Surv(AVAL, Status) ~ .") %>%
  forest_model()

forest_model(visR_HR)

surv_HR <- coxph(Surv(AVAL, Status) ~ ., lung_cohort)

summary(coxph(Surv(start, stop, event) ~ x, test2)) 

coxph(Surv(time, status) ~ ., data = pretty_lung)

my_func <- function(
  x,
  update_formula = NULL,
  ...
){
  
  # Update formula ----------------------------------------------------------
  
  if (!is.null(update_formula)){
    updated_object <- stats::update(x,  formula = eval(update_formula), evaluate = TRUE)
  } else updated_object <- x
  
  # Change Call -------------------------------------------------------------
  
  SurvCall <- as.list(updated_object$call)
  CoxArgs <- base::formals(survival::coxph)
  CoxCall <- append(quote(survival::coxph), SurvCall[names(SurvCall) %in% names(CoxArgs)])
  CoxCall <- append(CoxCall, list(...))
  
  # Tidy output -------------------------------------------------------------
  
  cox <- eval(as.call(CoxCall))
  
  return(cox)
}
