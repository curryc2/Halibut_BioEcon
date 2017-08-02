#' Read in and update halibut life history parameters form .xlsx. 
#'   Calculate size, mortality, and selectivity parameters at age, given sex and gear type. 
#'
#' @return halibut a list including all available growth, mortality, selectivity information for the simulated stock. 
#' @export
#'
read_update_params <- function() {
  
  load('data/halibut.rda')
  #=============================================================
  #Adjust halibut object values based on inputs from spreadsheet
  #Growth
  in.growth <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Growth')
  halibut$theta$linf <- as.numeric(in.growth[in.growth$Par=='linf',(2:3)])
  halibut$theta$vbk <- as.numeric(in.growth[in.growth$Par=='vbk',(2:3)])
  halibut$theta$to <- as.numeric(in.growth[in.growth$Par=='to',(2:3)])
  halibut$theta$a <- as.numeric(in.growth[in.growth$Par=='a',(2:3)])
  halibut$theta$b <- as.numeric(in.growth[in.growth$Par=='b',(2:3)])
  #Maturity
  in.maturity <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Maturity')
  halibut$theta$ahat <- as.numeric(in.maturity[in.maturity$Par=='ahat',(2:3)])
  halibut$theta$ghat <- as.numeric(in.maturity[in.maturity$Par=='ghat',(2:3)])
  #Mortality
  in.mortality <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Mortality')
  halibut$theta$m <- as.numeric(in.mortality[in.mortality$Par=='m',(2:3)])
  halibut$theta$A <- max(in.mortality[in.mortality$Par=='A',2:3])
  #Fishery Selectivity
  in.selex <- read.xlsx2('Halibut Model Inputs.xlsx', sheetName='FisherySelectivity', stringsAsFactors=FALSE, header=TRUE)
  halibut$MP$slx1 <- as.numeric(in.selex[in.selex$Par=='mu',c(2:5)])
  halibut$MP$slx2 <- as.numeric(in.selex[in.selex$Par=='sigma',c(2:5)])
  halibut$MP$slx3 <- as.numeric(in.selex[in.selex$Par=='gamma',c(2:5)])
  halibut$MP$slx4 <- as.numeric(in.selex[in.selex$Par=='plus.age',c(2:5)])  #Plus group age for selectivity
  halibut$MP$slim <- as.numeric(in.selex[in.selex$Par=='slim',c(2:5)])  #Minimum size for retention
  halibut$MP$ulim <- as.numeric(in.selex[in.selex$Par=='ulim',c(2:5)])  #Maximum size for retention
  halibut$MP$dmr <- as.numeric(in.selex[in.selex$Par=='dmr',c(2:5)])  #Discard Mortality Rate
  # halibut$MP$pscLimit <- as.numeric(in.selex[in.selex$Par=='pscLimit',c(2:5)])  #Prohibited species catch limit
  halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])
  halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])
  halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])
  
  #Recruitment Parameter Inputs
  in.rec <- read.xlsx2('Halibut Model Inputs.xlsx', sheetName='Recruitment', stringsAsFactors=FALSE)
  halibut$rec$steep <- as.numeric(in.rec$Value[in.rec$Par=='steep'])
  halibut$rec$sigma_rec <- as.numeric(in.rec$Value[in.rec$Par=='sigma_rec'])
  halibut$rec$ro  <- as.numeric(in.rec$Value[in.rec$Par=='ro'])
  
  #=========================================
  halibut <- getSelectivities(halibut)

  return(halibut)
}
