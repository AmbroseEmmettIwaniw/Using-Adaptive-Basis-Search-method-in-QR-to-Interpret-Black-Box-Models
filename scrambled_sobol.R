source("C:/Users/ambro/Documents/rsobol.R")

scrambled_sobol <- function(m=10,s=5,seed=20669099){
  if (s <= 50){
    rsobol(fn="C:/Users/ambro/Documents/fiftysobol.col", m=m, s=s,rand=TRUE, type="nestu",M=32,seed=seed)
  }
  else {
    rsobol(fn="C:/Users/ambro/Documents/sobol_Cs.col", m=m, s=s,rand=TRUE, type="nestu",M=32,seed=seed)
  }
}