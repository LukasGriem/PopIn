#-----------------------------------------------------------------------------------------#
################ PopIn - impact of linear Infrastructures on Populations ##################
#
# this file includes:
#   sigma2: dispersal variance
#   r1: population growth rate in favourable habitats
#   mlid: maximum linear infrastructure density above which a population cannot persist
#   mps: minimum size of a favourable habitat patch delimited by linear infrastructures 
#     below which a population cannot persist
#   indmodel: spatially-explicit individual-based model of population dynamics
#   agemodel: spatially-implicit stochastic age-based model stochastic of population dynamics
#
# function agemodel and function indmodel:
#    can receive an object (list) or a file with the values of the parameters
#    can run one or several replicates
#
#
#
#
#-----------------------------------------------------------------------------------------#
############################## sigma2 - dispersal variance ################################



# needs sigma_median, u

sigma2 <- function(sigma_median = 7*sqrt(129.89), u = 1/16) 
{ 
  ( (sigma_median)^2 / (1.18)^2 ) * ( 1 / (1/u) )
}

# sigma2: dispersal variance (km2/year)
# sigma_median: dispersal median (km/generation)
# u: mortality rate (per year)
# u can be computed from mean life span (mls) as: u = 1/mls
# sigma_median can be computed from home range (HR; km2) as: sigma_median = 7*sqrt(HR)
# references: Bowman et al. 2002, Pereira and Daily 2006, Ceia-Hasse et al. 2017
# default values for Puma concolor (Jones et al. 2009, Pereira and Daily 2006)
# sigma2() 
## [1] 285.6852
#
#################################### end of sigma2 ########################################
#-----------------------------------------------------------------------------------------#
#
#
#
#
#-----------------------------------------------------------------------------------------#
################### r1 - population growth rate in favourable habitats ####################

# needs b, B, s, maxage

r1 <- function(b = 1.75, B = 2.75, s = (1-(1/16)) , maxage = 80) 
{
  i <- seq(B, maxage)
  
  r <- uniroot(function(x) { 1-(sum(exp(-x*i)*(s^i)*b))  }, 
               interval = c(-10,10), 
               extendInt = "yes")
  
  return(r$root)
}

# b: 0.5*litter size
# B: breeding age
# s: survival rate
# maxage: maximum age
# see Pereira et al. 2004 for details
# default values for Puma concolor (Pereira and Daily 2006), where
#  s = 1 - u => s = 1 - (1/mls)
#  maxage = 5 * mls
#  with u = natural mortality rate and mls = mean life span
# r1() 
## [1] 0.4631507

###################################### end of r1 ##########################################
#-----------------------------------------------------------------------------------------#
#
#
#
#
#-----------------------------------------------------------------------------------------#
##################### mlid - maximum linear infrastructure density ########################

# needs r1, r0, li.width

mlid <- function(r1 = 0.4631507, r0 = -(1/16)*10^3, li.width = 0.01)
{
  Dmax <- r1/(r1 + abs(r0))/li.width

  cat("Dmax = ", Dmax, "km/km2","\n")
  
  return(c(Dmax))
}

# r1 and r0 values for Puma concolor

#################################### end of mlid ##########################################
#-----------------------------------------------------------------------------------------#
#
#
#
#
#-----------------------------------------------------------------------------------------#
############################### mps - minimum patch size ##################################

# needs r1, sigma2, alpha

mps <- function(r1 = 0.4631507, sigma2 = 285.6852, alpha = 1)
{
  Lmin <- pi *sqrt(sigma2/2/r1 * (1+1/alpha^2))
  Amin <- Lmin * Lmin * alpha
  
  cat("Lmin = ", Lmin, "km", "\n")
  cat("Amin = ", Amin, "km2", "\n")
  return(c(Lmin, Amin))
}

# r1 and sigma2 values for Puma concolor

#################################### end of mps ###########################################
#-----------------------------------------------------------------------------------------#
#
#
#
#
#-----------------------------------------------------------------------------------------#
####################################### indmodel ##########################################

# spatially-explicit individual-based population dynamics model - called from Cpp

#' @useDynLib PopIn
#' @importFrom Rcpp sourceCpp

indmodel <- function(land_r="", param="", file="")
{
  
  start.time <- date()
  
  if (isTRUE(identical(param, "")) && isTRUE(identical(file, "")))
  {
    # if the user does not provide a list or a file with the parameter values
    # use default file
    file.i <- system.file("extdata", "pim.txt", package="PopIn")
    con.i <- file(file.i, "r")
    while ( TRUE ) {
      line <- readLines(con.i, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      eval(parse(text = line))
    }
    close(con.i)
  }
  
  else if (isFALSE(identical(par, "")) && isFALSE(identical(file, "")) )
  { 
    # if the user provides a list AND a file with the parameter values
    # return warning message
    cat("too much information... use a list OR a file with the parameter values", "\n")
    return(NULL)
  }
  
  else if(any(param!="")){
    # if the user provides a list with the parameter values
    # read user list
    n_rep <- param$n_rep
    hab_aff <- param$hab_aff
    n_steps <- param$n_steps 
    init_population <- param$init_population 
    hr_size <- param$hr_size 
    birth_rate <- param$birth_rate
    breeding_age <- param$breeding_age
    survival <- param$survival
    distance_weight <- param$distance_weight 
    dispersal_distance <- param$dispersal_distance
    dispersal_mode <- param$dispersal_mode
    sink_avoidance <- param$sink_avoidance
    neigh_avoidance <- param$neigh_avoidance
    sink_mortality <- param$sink_mortality
    plotYN <- param$plotYN
    writeYN <- param$writeYN
    dir_name <- param$dir_name
    file_name <- param$file_name
    disturb_radius <- param$disturb_radius
    ext_point <- param$ext_point
  }
  
  else{
    # if the user provides a file with the parameter values
    # read user file
    con <- file(file, "r")
    while ( TRUE ) {
      line <- readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      eval(parse(text = line))
    }
    close(con)
  }
  
  land_r <- as.matrix(land_r)
  if(land_r[1,1]==""){   # matrix by default
    nrow=200; ncol=200
    land_r <- matrix(1,nrow,ncol)
    
####################################### indmodel ##########################################
    land_r[10,] <- 0; land_r[20,] <- 0; land_r[30,] <- 0; land_r[40,] <- 0
    land_r[50,] <- 0; land_r[60,] <- 0; land_r[70,] <- 0; land_r[80,] <- 0
    land_r[90,] <- 0; land_r[100,] <- 0; land_r[110,] <- 0; land_r[120,] <- 0; 
    land_r[130,] <- 0; land_r[140,] <- 0; land_r[150,] <- 0; land_r[160,] <- 0; 
    land_r[170,] <- 0; land_r[180,] <- 0; land_r[190,] <- 0

    land_r[,10] <- 0; land_r[,20] <- 0; land_r[,30] <- 0; land_r[,40] <- 0
    land_r[,50] <- 0; land_r[,60] <- 0; land_r[,70] <- 0; land_r[,80] <- 0
    land_r[,90] <- 0; land_r[,100] <- 0; land_r[,110] <- 0; land_r[,120] <- 0; 
    land_r[,130] <- 0; land_r[,140] <- 0; land_r[,150] <- 0; land_r[,160] <- 0; 
    land_r[,170] <- 0; land_r[,180] <- 0; land_r[,190] <- 0
  
    # replace the values of the matrix by the values of hab_aff
    land_r <- mapvalues(land_r, sort(unique(as.vector(land_r))), hab_aff[1:length(hab_aff)])
    
    # transpose matrix
    land_r <- t(land_r)
    
  } else {
    dimM <- dim(land_r)  # if matrix is given by user; determine dims and transpose matrix
    nrow <- dimM[1]
    ncol <- dimM[2]
    
    land_r <- mapvalues(land_r, sort(unique(as.vector(land_r))), hab_aff[1:length(hab_aff)])
    land_r <- t(land_r)
  }
  
  # run replicates and call the model from Cpp (indmodelseC)
  
  if(writeYN == T){
    system(paste("mkdir", dir_name))
  }
  else{}
  
  imse_rep <- NULL
  
  for(i in 1:n_rep){
    
    file_name_i <- paste(dir_name,"/",file_name,"_",i, sep="")
    
    imse_rep <- cbind(imse_rep, 
                      replicate(1, 
                                indmodelseC(land_r = land_r, nrow = nrow, ncol = ncol, n_steps = n_steps, 
                                            init_population = init_population, hr_size = hr_size, 
                                            birth_rate = birth_rate, breeding_age = breeding_age, 
                                            survival = survival, distance_weight = distance_weight, 
                                            dispersal_distance = dispersal_distance, dispersal_mode = dispersal_mode, 
                                            sink_avoidance = sink_avoidance, neigh_avoidance = neigh_avoidance, 
                                            sink_mortality = sink_mortality, file_name = file_name_i, 
                                            ext_point = ext_point, disturb_radius = disturb_radius)
                                )
                      )
  }  

  # imse_rep contains the population size at each time step, for each replicate

  # get population sizes from imse_rep
  popsize <- imse_rep[1,]
  
  ## time_steps <- imse_rep[2,]
  
  # create a matrix with the population size at each time step (columns), for each 
  #  replicate that is run (lines)
  popsize_matrix <-  matrix(unlist(popsize), ncol = n_steps+1, byrow = TRUE)
  
  # average population size in each time step
  if(n_rep>1){
    pop_size_av_time_step <- colMeans(popsize_matrix[,])
  }
  else{
    pop_size_av_time_step <- popsize_matrix[,]
  }
  
  # average population size at the end of the simulation
  pop_size_av_end <- mean(popsize_matrix[,n_steps+1])
  
  # probability of population extinction and probability of population persistence
  
  # number of populations (replicates) extinct at the last time step
  n_pop_ext <- sum(popsize_matrix[,n_steps+1]==0)
  
  # probability of population extinction
  prob_ext <- n_pop_ext / n_rep
  
  # probability of population persistence
  prob_persist <- 1 - prob_ext
  
  # average time to extinction
  time_ext_byrep <- apply(popsize_matrix, 1, function(n) sum(n!=0))
  
  time_ext_av <- mean(time_ext_byrep)
  
  # cat("start.time = ", start.time, "\n")
  # cat("end.time = ", date(), "\n", "\n")
  
  if(plotYN==TRUE){
    # plot population size as a function of time (steps)
    plot(seq(1,length(popsize_matrix[1,]),1), 
         pop_size_av_time_step, 
         xlim=c(0,length(popsize_matrix[1,])), 
         ylim=c(0,max(pop_size_av_time_step)), 
         type="l", xlab="Time steps", ylab="Population size")
  }
  
  # return results
  return(list(# "popsize" = popsize, 
              "popsize_av_end" = pop_size_av_end, 
              "prob_ext" = prob_ext, 
              "prob_persist" = prob_persist, 
              "time_ext_av" = time_ext_av, 
              "popsizematrix" = popsize_matrix,
              "popsize_av_time_step" = pop_size_av_time_step))
}

.onUnload <- function (libpath)
{
  library.dynam.unload("PopIn", libpath)
}
################################### end of indmodel #######################################
#-----------------------------------------------------------------------------------------#
#
#
#
#
#-----------------------------------------------------------------------------------------#
####################################### agemodel ##########################################

# spatially-implicit stochastic age-based model of population dynamics

agemodel <- function(n.rep, param="", file=""){
    
  start.time <- date()
    
  #######################################################################################
  #-------------------------------------------------------------------------------------#
  # define function agemodel_basic (it is called when n.rep = 1) 
  #  and function agemodel_rep (it is called when when n.rep > 1)
    
  # agemodel_basic: agemodel with n.rep = 1
    
  agemodel_basic <- function(years,transient1, transient2, transient3, 
                             f.mort, a, b, 
                             n.ages, min.age, repr.period, months.of.recruitment, 
                             prob.recr.month.m, prob.recr.month.sd, 
                             prop.recr.age.m, prop.recr.age.sd, 
                             clutch.size.m, clutch.size.sd, 
                             jsmm, jssdm, jsmf, jssdf,
                             ysmm, yssdm, ysmf, yssdf, 
                             asmm, assdm, asmf, assdf, 
                             plotTF, writeTF, dirname.basic)
  {
      
    # in case there is a problem the user may have feedback immediately
    #    currently not active (i.e., files are being overwritten)
    if(writeTF == T){
      dir <- paste("Res_",dirname.basic,sep="")
      if(file.exists(dir)==T){
        #      stop("\nRes' directory already exists!\nChoose a different 
        #           name for the results' directory.\n", call. = F)
      }
      else{
        dir.create(dir)
      }
    }
      
  #######################################################################################
  ##### auxiliary functions #############################################################
      
     est.annual.road.mortality <- function(N,N.year,f.mort){
        
        # distribute the dead among months proportionally 
        # to the abundance in the previous year:
        
        N.K <- rep(0,12)
        number.dead <-round(f.mort * N)
        prop <- N.year/sum(N.year)
        for(i in 1:12){
          n <- runif(number.dead)
          N.K[i] <- length(n[n<=prop[i]])
        }
        
        return(N.K)
      }
      
  #######################################################################################
      
      est.road.mortality <- function(N,N.mort){
        
        # in this case it distributes the number of dead individuals
        # by a matrix with 12 columns (months) and n.ages rows
        
        N.n <- N
        
        n <- runif(N.mort)
        for(k in 1:N.mort){
          
          # now randomly chooses the individuals
          csN <- cumsum( t(N/sum(N.n)) )
          
          pos <- min(which(csN >= n[k]))
          
          p.i <- ceiling(pos/12)  
          p.j <- (pos - 1)%%12 +1 
          
          N.n[p.i,p.j] <- N.n[p.i,p.j] - 1 
        }
        
        N <- N - N.n
        
        N[N<0] <- 0 # just in case 
        return(N)
      }
      
  #######################################################################################
      
      recr <- function(N, prob.recr.month.m, prob.recr.month.sd, prop.recr.age, 
                       a, b, clutch.size.m, clutch.size.sd){
        
        # Add all individuals in the same age class that contribute to recruitment
        
        n.recr <- apply(N,1,sum)
        
        # Convert the logit values into a linear probability
        pr <- rnorm(1,prob.recr.month.m, prob.recr.month.sd)
        prob.recr.month <- exp(pr)/(1+exp(pr))
        
        # Proportion breeding per age
        pba <- apply(prop.recr.age,1,function(x) rbeta(1,x[1],x[2]))
        n.recr <- pba * n.recr
        
        no.pairs <- trunc(sum(n.recr)/2 * prob.recr.month)
        clutch <- rnorm(no.pairs,clutch.size.m,clutch.size.sd)
        S <- sum(clutch[clutch > 0])
        
        no.recr <- a*S/(b+S)        # Beverton-Holt relationship
        
        return(no.recr)
      }
      
  #######################################################################################
      
      natural.mortality <- function(alp.m, bet.m, alp.f, bet.f){
        
        males   <- rbeta(1,alp.m,bet.m)
        females <- rbeta(1,alp.f,bet.f)
        
        males <- -log(males)
        females <- -log(females)
        
        mort <- (males + females)/2/12  # average and PER MONTH
        
        return(mort)
      }
      
  #######################################################################################
      
      beta.alpha.beta <- function(m,s){
        
        # given the mean (m) and the standard deviation (s) this function calculates
        # the alpha and beta parameters of the beta distribution
        
        common <- (m-m^2)/s^2 - 1
        alp <- m*common
        bet <- (1-m)*common
        
        alp[m==1] <- 1
        bet[m==1] <- 0
        
        return(cbind(alp,bet))
      }
      
  #######################################################################################
      
    translate.month <- function(month.recr){
      
      n <- length(month.recr)
        
      if(month.recr[1] == "all" | month.recr[1] == "All"){
        months <- rep(1,12) 
      }
      else{
        months <- rep(0,12)
        for(i in 1:n){
         if(month.recr[i] == "January"   | month.recr[i] == "january")   months[1] <- 1
         if(month.recr[i] == "February"  | month.recr[i] == "february")  months[2] <- 1
         if(month.recr[i] == "March"     | month.recr[i] == "march")     months[3] <- 1
         if(month.recr[i] == "April"     | month.recr[i] == "april")     months[4] <- 1
         if(month.recr[i] == "May"       | month.recr[i] == "may")       months[5] <- 1
         if(month.recr[i] == "June"      | month.recr[i] == "june")      months[6] <- 1
         if(month.recr[i] == "July"      | month.recr[i] == "july")      months[7] <- 1
         if(month.recr[i] == "August"    | month.recr[i] == "august")    months[8] <- 1
         if(month.recr[i] == "September" | month.recr[i] == "september") months[9] <- 1
         if(month.recr[i] == "October"   | month.recr[i] == "october")   months[10] <- 1
         if(month.recr[i] == "November"  | month.recr[i] == "november")  months[11] <- 1
         if(month.recr[i] == "December"  | month.recr[i] == "december")  months[12] <- 1
        }
      }
      return(months)
    }
      
  #######################################################################################
      
      is.it.viable <- function(months.recr,min.age,repr.period,n.ages){
        
        # what we are checking here is that if there are any months in common:
        # 1st month of repr + min.age + n* repr.period <=> months.recr
        
        months.recr <- which(months.recr == 1)
        m <- floor(min.age/12)
        months <- months.recr + m*12
        months <- c(months, months + 12)
        
        first.age <- months.recr[1] + min.age -1   
        ages   <- seq(first.age, first.age+12,repr.period)
        
        res <- length(intersect(months,ages))
        
        return(res)
      }
      
  ##### end of auxiliary functions ######################################################
  #######################################################################################
  #
  # "years" and "transient1, 2, 3" are set in years, however,
  # the time unit used in the code is the month.
  #
  # - the first transient (T1) ensures that the age classes stabilize;
  # - the second transient (T2) is used to estimate the average population size
  # without road mortality;
  # - the third transient (T3) ensures that the population "stabilizes" after
  # road mortality starts being applied. Schematically:
  #
  #         ^<-T1-><-   T2   -><-T3->
  #         |
  #   pop.  | *         
  #   size  |*    ************* 
  #         |  * *             *
  #         |   *               *
  #         |                    *    
  #         |                     ****************
  #         |                          
  #      ----------------------------------------------->
  #         0                time                T
  #
  #
  # PDFs with values between 0 and 1 are modeled with beta distributions. 
  # The function beta.alpha.beta transforms the mean and sd into the 
  # corresponding alpha and beta parameters of the beta distribution.
  #
  # - Mortality has two components: natural mortality and road mortality.
  # Notice that the function natural.mortality already divides the sampled value
  # by 12, therefore taking into account that it is per month
  #
  # f.mort <- is the fraction of animals killed per year.
  #
  # RECRUITMENT. There are three aspects to be considered: 
  #
  #  - min.age: age at maturity, that is, the minimum age at which an animal can 
  #    reproduce. It is given in months, so if it is =13 it means it starts
  #    reproducing in the first month of its second year.
  #  - repr.period: how often a female reproduces; interval between clutches/litters (in months).
  #    If =5 it means it reproduces every 5 months.
  #    If =1 it means it reproduces every month.
  #  - months.of.recruitment: the months where animals reproduce.
  #    if = c("June","August") or = c("june","august") then it can reproduce only 
  #    in June and August. 
  #
  #  - In addition there are two other aspects to be taken into account:
  #  
  #    - In those months where the population can reproduce, which proportion
  #      reproduces: variables (vectors) prob.recr.month.m(sd) (in a logit scale)
  # 
  #    - What is the proportion per age that reproduces: variables (vectors)
  #      prop.recr.age.m(sd) (in a linear scale)
  #
    
  # recalculate the transients:
      
      transient2 <- transient1 + transient2
      transient3 <- transient2 + transient3
      
  #######################################################################################
      # Initialize the population:
      
      # n.ages is the number of age classes, ie, the maximum number of
      # years an individual can live.
      
      # N is a matrix containing the number of individuals. It has 12 columns and
      # n.ages rows. So the age of the individuals is identified per month (and 
      # not simply per year). For instance, row 2 and column 3 contains the
      # number of individuals that are 1 year and 2 months old.
      
      N.ini <- 30           # N.ini is the initial population per age class
      # 30 is just an initial value
      N <- matrix(N.ini,nrow=n.ages,ncol=12)
      
      
  #######################################################################################
    # RECRUITMENT:
      
    # min.age 
    # min.age is the age of first breeding; only animals this 
    # age or older (IN MONTHS) can reproduce. 
      
    # repr.period
    # repr.period is the reproduction interval for a female
    # that is, the interval in months where a female gives birth
      
    # months.of.recruitment
    # These months are translated into 0s (no recruitment) and 1s (recruitment):
    months.recr <- translate.month(months.of.recruitment)
      
    # proportion breeding per MONTH (mean and stdev):
    # prob.recr.month.m 
    # prob.recr.month.sd
    # if only a fraction of the mature reproduce in a specific MONTH change
    # the vectors of the mean and sd 
    # (Notice that there is a certain redundancy in the number of variables.
    # We could be using prob.recr.month instead of months.recr.)
      
    # Values correspond to logit, hence the 100 for the mean => prob ~ 1.
      
    # proportion breeding per AGE (mean and stdev):
    # (notice that prob.recr.month is per month while prop.recr.age is per age)  
    # prop.recr.age.m
    # prop.recr.age.sd
    prop.recr.age <- beta.alpha.beta(prop.recr.age.m,prop.recr.age.sd)  # cbind(alpha,beta)
      
    # clutch size (mean and stdev) # assuming normal distribution!
    # clutch.size.m
    # clutch.size.sd
      
    # parameters of the Beverton-Holt relationship:
    #  a
    #  b
    # recall that the number of fledgings has to be smaller than the number of
    #  eggs, that is why Rickers model does not work. For the Beverton-Holt
    #  model to be reasonable: a <= b
    # see Borda-de-Água et al. 2014 for more details
      
  #######################################################################################
      # SURVIVAL RATES: 
      # The mean and stdev are transformed into the alpha and beta parameters
      # of the beta distribution (function beta.alpha.beta)
      
      # juvenile survival rate
      res <- beta.alpha.beta(jsmm,jssdm)      # MALES
      # jsmm: mean of juvenile survival rate for males
      # jssdm: stdev of juvenile survival rate for males
      jnmm.alp <- res[1]; jnmm.bet <- res[2]
      
      res <- beta.alpha.beta(jsmf,jssdf)      # FEMALES
      # jsmf: mean of juvenile survival rate for females
      # jssdf: stdev of juvenile survival rate for females
      jnmf.alp <- res[1]; jnmf.bet <- res[2]
      
      # yearling survival rate 
      res <- beta.alpha.beta(ysmm,yssdm)     # MALES
      # ysmm: mean of yearling survival rate for males
      # yssdm: stdev of yearling survival rate for males
      yeam.alp <- res[1]; yeam.bet <- res[2]
      
      res <- beta.alpha.beta(ysmf,yssdf)     # FEMALES
      # ysmf: mean of yearling survival rate for females
      # yssdf: stdev of yearling survival rate for females
      yeaf.alp <- res[1]; yeaf.bet <- res[2]
      
      # adult survival rate
      res <- beta.alpha.beta(asmm,assdm)     # MALES
      # asmm: mean of adult survival rate for males
      # assdm: stdev of adult survival rate for males  
      adum.alp <- res[1]; adum.bet <- res[2]
      
      res <- beta.alpha.beta(asmf,assdf)     # FEMALES
      # asmf: mean of adult survival rate for females
      # assdf: stdev of adult survival rate for females   
      aduf.alp <- res[1]; aduf.bet <- res[2]
      
      
      # End of inputing parameters
      
  #######################################################################################
      
      # Initialize matrices and vectors:
      
      N.total  <- rep(0,12*years)  # total population size per month/iteration
      N.annual.jan <- rep(0,years) # population size in January
      N.annual.dec <- rep(0,years) # population size in December
      N.killed <- rep(0,years)
      road.mortality <- matrix(0,nrow=n.ages,ncol=12)
      
      f.mort.12 <- -1                                                
      
      # Identify the months that contribute to the recruitment,
      #  during the entire life span of an individual (r.m), 
      #  without taking into account which are the yearly months of recruitment
      # (These two sets of information will be combined later on in the program)
      
      r.m <- seq(min.age,12*n.ages,repr.period)
      #print("r.m")
      #print(r.m)
      
      # Based on this information we construct a matrix of 0s (no recruitment)
      #  and 1s (recruitment).
      # This matrix (mat.recr) has 12 (months) columns and n.ages (years) rows.
      # In the months where recruitment is allowed this matrix "multiplies" N,
      # thus selecting those individuals that can contribute to reproduction.
      
      sel.row <- ceiling(r.m/12)
      sel.col <- (r.m-1) %% 12 + 1
      age.recr <- unique(sel.row)
      
      mat.recr <- matrix(0,ncol=12,nrow=n.ages)
      for(i in age.recr)
        mat.recr[i,sel.col[sel.row==i]] <- 1
      
      #print("mat.recr")
      #print(mat.recr)
      can.repr <- mat.recr 
      
      #print("can.repr")
      #print(can.repr)
      
      # THE MOST IMPORTANT CYCLES:
      # for a given number of years, advance the population monthly
      
      for(y in 1:years){   # Years
        
        for(m in 1:12){   # Months
          
          #cat("y= ", y," m= ",m,"\n")
          #print("BEFORE TIME ADVANCEMENT")
          #print("N")
          #print(N)
          
        # Estimate mortality at different age classes; notice this is a random process.
          
          juv.nat.mort <- natural.mortality(jnmm.alp,jnmm.bet,jnmf.alp,jnmf.bet)
          yea.nat.mort <- natural.mortality(yeam.alp,yeam.bet,yeaf.alp,yeaf.bet)
          adu.nat.mort <- natural.mortality(adum.alp,adum.bet,aduf.alp,aduf.bet)
          
          Z <- c(juv.nat.mort,yea.nat.mort,rep(adu.nat.mort,(n.ages-2)))
          
          # TIME ADVANCEMENT
          # Advances one month the entire population.
          
          # first the ones that were in the (year + 12 month) become 1 year old
          # (we keep this information); the oldest ones die:
          one.year.older <- N[1:(n.ages-1),12]
          
          # all the other ones advance one month
          N[,2:12] <- N[,1:11]
          
          # now update the ones that become one year old
          N[2:n.ages,1] <- one.year.older
          
          N[1,1] <- 0 
          
          #print("TIME ADVANCEMENT")
          #print("N")
          #print(N)
          
          # MORTALITY
          # Apply natural mortality
          # Estimate mortality at different age classes; notice this is a random process.
          
          juv.nat.mort <- natural.mortality(jnmm.alp,jnmm.bet,jnmf.alp,jnmf.bet)
          yea.nat.mort <- natural.mortality(yeam.alp,yeam.bet,yeaf.alp,yeaf.bet)
          adu.nat.mort <- natural.mortality(adum.alp,adum.bet,aduf.alp,aduf.bet)
          
          Z <- c(juv.nat.mort,yea.nat.mort,rep(adu.nat.mort,(n.ages-2)))
          N <- N * exp(-Z)
          
          N <- round(N); N[N<0] <- 0; nt <- sum(N); if(nt <= 0) {break}
          
          # Apply road mortality
          # Road mortality only occurs after transient 2 (defined in years)
          
          if(y > transient2 & road.mortality[m] > 0) {
            
            if(road.mortality[m] > nt) N.mort <- nt
            else N.mort <- road.mortality[m] 
            
            N.k <- est.road.mortality(N,N.mort)
            N <- N-N.k
            
            N.killed[y] <- N.killed[y] + sum(N.k)
            
            N <- round(N); N[N<0] <- 0; nt <- sum(N); if(nt <= 0) {break}
          }
          if(nt <= 0) {break}
          
          #print("AFTER MORTALITY")
          #print("N")
          #print(N)
          
          # RECRUITMENT 
          # advances the can.repr matrix (but only it it exists)
          if(months.of.recruitment[1] != "all"){
            if(length(can.repr) != 1){
              one.year.older <- can.repr[1:(n.ages-1),12]
              can.repr[,2:12] <- can.repr[,1:11]
              can.repr[2:n.ages,1] <- one.year.older
              can.repr[1,1] <- 0
            }
          }
          # identifies the age classes that can reproduce
          can.repr <- can.repr + mat.recr
          
          # Introduce recruits into the population. 
          # Only the ones that survive this month contribute to the recruitment. 
          # But only in those months where recruitment is allowed (months.recr==1). 
          # Notice the "N*can.repr", this selects the individuals that can 
          # contribute to recruitment
          
          if(months.recr[m] == 1){
            #print("N*can.repr")
            #print(N*can.repr)
            N[1,1] <- recr(N*can.repr, prob.recr.month.m[m], prob.recr.month.sd[m], 
                           prop.recr.age, a, b, clutch.size.m, clutch.size.sd)
            N[1,1] <- round(N[1,1]); nt <- sum(N); if(nt <= 0) {break}
            
            # it assumes that all the females that could reproduce, reproduce
            can.repr <- 0
            #print("AFTER RECRUITMENT")
            #print("N")
            #print(N)
            #print("can.repr")
            #print(can.repr)
          }
          
          # now collect information
          
          N.total[(y-1)*12+m] <- nt
          
          if(m==12 & y>= (transient2-1)){
            
            f.mort.12 <- c(f.mort.12,1-N.total[(y-1)*12+m]/N.total[(y-2)*12+m])
            
            # Estimate the mortality to be applied to the next year
            road.mortality <- est.annual.road.mortality(nt,N.total[((y-1)*12+1):((y-1)*12+12)], f.mort)
          }
          
        } # end monthly cycle
        
        # collect information to calculate the average before road mortality starts
        
        N.annual.jan[y] <- N.total[(y-1)*12 + 1]      # January
        N.annual.dec[y] <- N.total[(y-1)*12 + 12]     # December
        
        if(nt <= 0) {break}
        
      } # end annual cycle
      
      no.it <- y # number of iterations
      
  #######################################################################################
  ##### PLOTS ###########################################################################
      
      if(plotTF== T){
        
        dev.new() # quartz()
        
        layout(matrix(c(1,2),nrow=2))
        
        plot(N.annual.jan[1:no.it],type="l",
             ylim=c(0,max(N.annual.jan)),
             ylab = "Population size January", xlab = "years")
        leg.txt <- c(paste("mean pop. size=",round(mean(N.annual.jan[transient3:no.it]),2)))
        legend("topright",leg.txt,lty=c(1),col=c("black"), bty="n")
        
        plot(N.annual.dec[1:no.it],type="l", 
             ylim=c(0,max(N.annual.dec)),
             ylab = "Population size December", xlab = "years")
        leg.txt <- c(paste("mean pop. size=",round(mean(N.annual.dec[transient3:no.it]),2)))
        legend("topright",leg.txt,lty=c(1),col=c("black"), bty="n")
      }
      
      if(no.it > transient3 & plotTF == T){
        
        dev.new() # quartz()
        
        plot(transient3:no.it,N.killed[transient3:no.it], 
             ylim=c(0,max(N.killed[transient3:no.it])), type="l",
             ylab = "Number killed", xlab = "years", main="")
        lines(c(transient3,transient3),c(0,max(N.killed[transient3:no.it])),col="red", lty=2)
        leg.txt <- c(paste("mean killed=",round(mean(N.killed[transient3:no.it]),3))) 
        legend("topright",leg.txt,lty=c(1),col=c("black"), bty="n")
      }
      
  ##### Histograms ######################################################################
      
      if(plotTF== T){
        dev.new(width=10,height=5) # quartz(width=10,height=5)
        layout(matrix(1:2,ncol=2))
        
        N.kil <- N.killed[N.killed != 0]
        
        res <- hist(log(N.kil),nclass=20,plot=T, xlab="log(Number Killed)", probability=T, main="")
        
        mids   <- res$mids[res$density != 0]
        density <- log(res$density[res$density!= 0])
        
        plot(mids,density,xlab="log(Number Killed)",ylab="log(Density)")
        
        div <- 5
        densityR <- density[div:length(mids)]
        midsR   <- mids[div:length(mids)]
        
        res <- lm(densityR ~ midsR)
        interc <- res$coefficients[[1]]
        slope  <- res$coefficients[[2]]
        lines(midsR, interc+slope*midsR, col="red")
        
        txt1 <- paste("slope =", round(slope,2))
        txt2 <- paste("intercept =", round(interc,2))
        
        densityL <- density[1:div]
        midsL   <- mids[1:div]
        
        res <- lm(densityL ~ midsL)
        interc <- res$coefficients[[1]]
        slope  <- res$coefficients[[2]]
        
        legend("topright",c(txt1,txt2),bty="n")
        
      } 
  ##### End of Histograms ###############################################################
      
      f.mort.12 <- f.mort.12[f.mort.12 >= 0]
      N.annual.jan <- N.annual.jan[N.annual.jan != 0]
      N.annual.dec <- N.annual.dec[N.annual.dec != 0]
      
      if(length(N.annual.jan)==0) N.annual.jan <- rep(0,years) 
      if(length(N.annual.dec)==0) N.annual.dec <- rep(0,years)
      
      res1  <- no.it
      res2  <- median(N.annual.jan[(transient1+1):transient2])      # no road mortality
      res8  <-     sd(N.annual.jan[(transient1+1):transient2])      # no road mortality
      res13 <-   mean(N.annual.jan[(transient1+1):transient2])      # no road mortality
      res18 <-    max(N.annual.jan[(transient1+1):transient2])      # no road mortality
      res19 <-    min(N.annual.jan[(transient1+1):transient2])      # no road mortality
      
      res3  <- median(N.annual.dec[(transient1+1):transient2])      # no road mortality
      res9  <-     sd(N.annual.dec[(transient1+1):transient2])      # no road mortality
      res14 <-   mean(N.annual.dec[(transient1+1):transient2])      # no road mortality
      res20 <-    max(N.annual.dec[(transient1+1):transient2])      # no road mortality
      res21 <-    min(N.annual.dec[(transient1+1):transient2])      # no road mortality
      
      if(y > transient3){ res4  <- median(N.annual.jan[(transient3+1):no.it]);
      res10 <-     sd(N.annual.jan[(transient3+1):no.it]);
      res15 <-   mean(N.annual.jan[(transient3+1):no.it]);
      res22 <-    max(N.annual.jan[(transient3+1):no.it]);
      res23 <-    min(N.annual.jan[(transient3+1):no.it])}
      else              { res4  <- -1; res10 <- -1; res15 <- -1; res22 <- -1; res23 <- -1;}
      
      if(y > transient3){ len.vec <- length(N.annual.dec);
      res5  <- median(N.annual.dec[(transient3+1):len.vec]);
      res11 <-     sd(N.annual.dec[(transient3+1):len.vec]);
      res16 <-   mean(N.annual.dec[(transient3+1):len.vec]);
      res24 <-    max(N.annual.dec[(transient3+1):len.vec]);
      res25 <-    min(N.annual.dec[(transient3+1):len.vec])}
      else              { res5  <- -1; res11 <- -1; res16 <- -1; res24 <- -1; res25 <- -1;}
      
      if(y > transient3){ res6  <- median(N.killed[(transient3+1):no.it]);
      res12 <-     sd(N.killed[(transient3+1):no.it]);
      res17 <-   mean(N.killed[(transient3+1):no.it]);
      res26 <-    max(N.killed[(transient3+1):no.it]);
      res27 <-    min(N.killed[(transient3+1):no.it])}
      else              { res6  <- -1; res12 <- -1; res17 <- -1; res26 <- -1; res27 <- -1;}
      res7 <- mean(f.mort.12)
      
      if(plotTF==F){
        return(c(res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, 
                 res11, res12, res13, res14, res15, res16, res17, res18, res19, res20, 
                 res21, res22, res23, res24, res25, res26, res27))
      }
      
      #   return(c(mean(N.annual.dec[transient3:no.it]), mean(N.annual.dec[transient3:no.it]), 
      #            mean(N.killed[transient3:no.it]), median(N.annual.dec[transient3:no.it]), 
      #            median(N.annual.jan[transient3:no.it]), median(N.killed[transient3:no.it])))
      
  ##### Write files #####################################################################
      
      mjan <- matrix(0,no.it,2)
      mjan[,1] <- seq(1,no.it,1)
      mjan[,2] <- N.annual.jan[1:no.it]
      
      mdec <- matrix(0,no.it,2)
      mdec[,1] <- seq(1,no.it,1)
      mdec[,2] <- N.annual.dec[1:no.it]
      
      mkill <- matrix(0,length(transient3:no.it),2)
      mkill[,1] <- transient3:no.it
      mkill[,2] <- N.killed[transient3:no.it]
      
      if(writeTF == T){
        write(t(mjan),paste(dir,"/pop_size_jan",sep=""), ncolumns=2) 
        write(t(mdec),paste(dir,"/pop_size_dec",sep=""), ncolumns=2) 
        write(t(mkill),paste(dir,"/n_killed",sep=""), ncolumns=2) 
      }
      
      cat("start.time = ", start.time, "\n")
      cat("end.time = ", date(), "\n", "\n")
      
  ############################################################################################
  #------------------------------------------------------------------------------------------#
  return(list("n.rep" = n.rep,
              "years" = years, 
              "transient1" = transient1, "transient2" = transient2, "transient3" = transient3, 
              "f.mort" = f.mort, "a" = a, "b" = b, 
              "n.ages" = n.ages, "min.age" = min.age, "repr.period" = repr.period,
              "months.of.recruitment" = months.of.recruitment, 
              "prob.recr.month.m" = prob.recr.month.m, 
              "prob.recr.month.sd" = prob.recr.month.sd,  
              "prop.recr.age.m" = prop.recr.age.m, "prop.recr.age.sd" = prop.recr.age.sd, 
              "clutch.size.m" = clutch.size.m, "clutch.size.sd" = clutch.size.sd, 
              "jsmm" = jsmm, "jssdm" = jssdm, "jsmf" = jsmf, "jssdf" = jssdf, 
              "ysmm" = ysmm, "yssdm" = yssdm, "ysmf" = ysmf, "yssdf" = yssdf, 
              "asmm" = asmm, "assdm" = assdm, "asmf" = asmf, "assdf" = assdf,
              "plotTF" = plotTF, "writeTF" = writeTF, "dirname" = dirname.basic))   
  #-----------------------------------------------------------------------------------------#
  ###########################################################################################
      
  }
  # end of agemodel_basic
  #-------------------------------------------------------------------------------------#
  #######################################################################################
  #
  #
  #######################################################################################
  #-------------------------------------------------------------------------------------#
  # agemodel_rep: agemodel with n.rep > 1
    
  agemodel_rep <- function(n.rep, years, transient1, transient2, transient3, 
                           a.val, f.values, 
                           n.ages, min.age, repr.period, months.of.recruitment, 
                           prob.recr.month.m, prob.recr.month.sd, 
                           prop.recr.age.m, prop.recr.age.sd, 
                           clutch.size.m, clutch.size.sd, 
                           jsmm, jssdm, jsmf, jssdf, 
                           ysmm, yssdm, ysmf, yssdf, 
                           asmm, assdm, asmf, assdf, 
                           plotTF, writeTF, dirname.rep)
    {
      
      # in case there is a problem the user may have feedback immediately
      #    currently not active (i.e., files are being overwritten)
      if(writeTF == T){
        dir <- paste("Res_",dirname.rep,sep="")
        if(file.exists(dir)==T){
          #      stop("\nRes' directory already exists!\nChoose a different 
          #           name for the results' directory.\n", call. = F)
        }
        else{
          dir.create(dir)
        }
      }
      
      t <- rep(0,n.rep)
      
      # transient3 <- 50  # if the pop. survives for more than 50 years collect info
      
      # recalculate the transients:
      
      transient1.n <- transient1
      transient2.n <- transient1.n + transient2
      transient3.n <- transient2.n + transient3
      
      # a.val
      # "a" and "b" are parameters of the Beverton-Holt relationship
      # "a.val" corresponds to the values that "a" will take
      # b will be set below to a.val[a] (condition: b >= a)
      # see Borda-de-Água et al. 2014 for more details
      
      n.a <- length(a.val)
      
      N.annual.jan  <- rep(0,n.a) # mean size at the end of the month before the 1st recruitment
      N.annual.jan.sd  <- rep(0,n.a)
      N.annual.jan.mean  <- rep(0,n.a)
      N.annual.jan.max   <- rep(0,n.a)
      N.annual.jan.min   <- rep(0,n.a)
      N.annual.jan.max.max <- rep(0,n.a)
      N.annual.jan.min.min <- rep(1e7,n.a)
      
      N.annual.dec <- rep(0,n.a)
      N.annual.dec.sd <- rep(0,n.a)
      N.annual.dec.mean <- rep(0,n.a)
      N.annual.dec.max  <- rep(0,n.a)
      N.annual.dec.min  <- rep(0,n.a)
      N.annual.dec.max.max <- rep(0,n.a)
      N.annual.dec.min.min <- rep(1e7,n.a)
      
      # f.values
      # f.values corresponds to the values that the fraction of
      # animals killed per year will take.
      
      n.f <- length(f.values)
      
      mean.t   <- matrix(0,nrow=n.a,ncol=n.f)
      median.t <- matrix(0,nrow=n.a,ncol=n.f)
      sd.t     <- matrix(0,nrow=n.a,ncol=n.f)
      max.t    <- matrix(0,nrow=n.a,ncol=n.f)
      min.t    <- matrix(0,nrow=n.a,ncol=n.f)
      prob.ext  <- matrix(0,nrow=n.a,ncol=n.f)
      f.mort.12 <- matrix(0,nrow=n.a,ncol=n.f)
      
      N.annual.jan.m      <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.sd   <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.mean <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.max <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.min <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.max.max <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.jan.m.min.min <- matrix(1e7,nrow=n.a,ncol=n.f)
      
      N.annual.dec.m      <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.sd   <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.mean <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.max  <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.min  <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.max.max  <- matrix(0,nrow=n.a,ncol=n.f)
      N.annual.dec.m.min.min  <- matrix(1e7,nrow=n.a,ncol=n.f)
      
      N.killed      <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.sd   <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.mean <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.max <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.min <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.max.max <- matrix(0,nrow=n.a,ncol=n.f)
      N.killed.min.min <- matrix(1e7,nrow=n.a,ncol=n.f)
      
      for (a in 1:n.a){
        for (f in 1:n.f){
          
          b.val <- a.val
          n.b <- length(b.val)
          b <- a.val[a]
          cat("a= ", a.val[a], "b= ", b, " f= ", f.values[f],"\n") 
          n.rep.m <- 0
          for (i in 1:n.rep){
            
            # call agemodel_basic()
            res <- agemodel_basic(years,transient1, transient2, transient3, 
                                  f.values[f], a.val[a], b, 
                                  n.ages, min.age, repr.period, months.of.recruitment, 
                                  prob.recr.month.m, prob.recr.month.sd, 
                                  prop.recr.age.m, prop.recr.age.sd, 
                                  clutch.size.m, clutch.size.sd, 
                                  jsmm, jssdm, jsmf, jssdf, 
                                  ysmm, yssdm, ysmf, yssdf, 
                                  asmm, assdm, asmf, assdf,
                                  plotTF=F, writeTF=F)
            
            t[i] <- res[1]
            
            # results for i > transient1 and i <= transient2
            
            N.annual.jan[a]      <- N.annual.jan[a] + res[2]
            N.annual.jan.sd[a]   <- N.annual.jan.sd[a] + res[8]
            N.annual.jan.mean[a] <- N.annual.jan.mean[a] + res[13]
            N.annual.jan.max[a]  <- N.annual.jan.max[a] + res[18]
            N.annual.jan.min[a]  <- N.annual.jan.min[a] + res[19]
            if(is.na(res[18]) == T) res[18] <- 0
            if(is.na(res[19]) == T) res[19] <- 1e7
            if(res[18] > N.annual.jan.max.max[a]) N.annual.jan.max.max[a] <- res[18]
            if(res[19] < N.annual.jan.min.min[a]) N.annual.jan.min.min[a] <- res[19]
            
            N.annual.dec[a]      <- N.annual.dec[a] + res[3]
            N.annual.dec.sd[a]   <- N.annual.dec.sd[a] + res[9]
            N.annual.dec.mean[a] <- N.annual.dec.mean[a] + res[14]
            N.annual.dec.max[a]  <- N.annual.dec.max[a] + res[20]
            N.annual.dec.min[a]  <- N.annual.dec.min[a] + res[21]
            if(is.na(res[20]) == T) res[20] <- 0
            if(is.na(res[21]) == T) res[21] <- 1e7
            if(res[20] > N.annual.dec.max.max[a]) N.annual.dec.max.max[a] <- res[20]
            if(res[21] < N.annual.dec.min.min[a]) N.annual.dec.min.min[a] <- res[21]
            
            if(res[1] == years){     # meaning the simulation reached the end
              
              # results after transient3 until years
              
              n.rep.m <- n.rep.m + 1  
              
              N.annual.jan.m[a,f]      <- N.annual.jan.m[a,f] + res[4]
              N.annual.jan.m.sd[a,f]   <- N.annual.jan.m.sd[a,f] + res[10]
              N.annual.jan.m.mean[a,f] <- N.annual.jan.m.mean[a,f] + res[15]
              N.annual.jan.m.max[a,f]  <- N.annual.jan.m.max[a,f] + res[22]
              N.annual.jan.m.min[a,f]  <- N.annual.jan.m.min[a,f] + res[23]
              if(is.na(res[22]) == T) res[22] <- 0
              if(is.na(res[23]) == T) res[23] <- 1e7
              if(res[22] > N.annual.jan.m.max.max[a,f]) N.annual.jan.m.max.max[a,f] <- res[22]
              if(res[23] < N.annual.jan.m.min.min[a,f]) N.annual.jan.m.min.min[a,f] <- res[23]
              
              N.annual.dec.m[a,f]      <- N.annual.dec.m[a,f] + res[5]
              N.annual.dec.m.sd[a,f]   <- N.annual.dec.m.sd[a,f] + res[11]
              N.annual.dec.m.mean[a,f] <- N.annual.dec.m.mean[a,f] + res[16]
              N.annual.dec.m.max[a,f]  <- N.annual.dec.m.max[a,f] + res[24]
              N.annual.dec.m.min[a,f]  <- N.annual.dec.m.min[a,f] + res[25]
              if(is.na(res[24]) == T) res[24] <- 0 
              if(is.na(res[25]) == T) res[25] <- 1e7 
              if(res[24] > N.annual.dec.m.max.max[a,f]) N.annual.dec.m.max.max[a,f] <- res[24]
              if(res[25] < N.annual.dec.m.min.min[a,f]) N.annual.dec.m.min.min[a,f] <- res[25]
              
              N.killed[a,f]      <- N.killed[a,f] + res[6]
              N.killed.sd[a,f]   <- N.killed.sd[a,f] + res[12]
              N.killed.mean[a,f] <- N.killed.mean[a,f] + res[17]
              N.killed.max[a,f]  <- N.killed.max[a,f] + res[26]
              N.killed.min[a,f]  <- N.killed.min[a,f] + res[27]
              if(is.na(res[26]) == T) res[26] <- 0
              if(is.na(res[27]) == T) res[27] <- 1e7
              if(res[26] > N.killed.max.max[a,f]) N.killed.max.max[a,f] <- res[26]
              if(res[27] < N.killed.min.min[a,f]) N.killed.min.min[a,f] <- res[27]
              
              f.mort.12[a,f] <- f.mort.12[a,f] + res[7]
            }
          }
          
          survived <- length(t[t > transient2.n])
          t <- t[t <= transient3.n & t > transient2.n]
          
          if(length(t) != 0){  # results for the populations that went extinct
            
            t <- t-transient2.n
            prob.ext[a,f] <- length(t)/survived
            mean.t[a,f] <- mean(t)
            median.t[a,f] <- median(t)
            sd.t[a,f]   <- sd(t)
            max.t[a,f]  <- max(t)
            min.t[a,f]  <- min(t)
          }
          
          if(n.rep.m > 0){     # results for the populations that survived to the end
            
            N.annual.jan.m[a,f]      <- N.annual.jan.m[a,f]/n.rep.m
            N.annual.jan.m.sd[a,f]   <- N.annual.jan.m.sd[a,f]/n.rep.m
            N.annual.jan.m.mean[a,f] <- N.annual.jan.m.mean[a,f]/n.rep.m
            N.annual.jan.m.max[a,f]  <- N.annual.jan.m.max[a,f]/n.rep.m
            N.annual.jan.m.min[a,f]  <- N.annual.jan.m.min[a,f]/n.rep.m
            
            N.annual.dec.m[a,f]      <- N.annual.dec.m[a,f]/n.rep.m
            N.annual.dec.m.sd[a,f]   <- N.annual.dec.m.sd[a,f]/n.rep.m
            N.annual.dec.m.mean[a,f] <- N.annual.dec.m.mean[a,f]/n.rep.m
            N.annual.dec.m.max[a,f]  <- N.annual.dec.m.max[a,f]/n.rep.m
            N.annual.dec.m.min[a,f]  <- N.annual.dec.m.min[a,f]/n.rep.m
            
            N.killed[a,f]      <- N.killed[a,f]/n.rep.m
            N.killed.sd[a,f]   <- N.killed.sd[a,f]/n.rep.m
            N.killed.mean[a,f] <- N.killed.mean[a,f]/n.rep.m
            N.killed.max[a,f]  <- N.killed.max[a,f]/n.rep.m
            N.killed.min[a,f]  <- N.killed.min[a,f]/n.rep.m
            
            f.mort.12[a,f] <- f.mort.12[a,f]/n.rep.m
          }
        }
      }
      
      # averaged values before road mortality
      
      N.annual.jan      <- N.annual.jan/(n.rep*n.f)
      N.annual.jan.sd   <- N.annual.jan.sd/(n.rep*n.f)
      N.annual.jan.mean <- N.annual.jan.mean/(n.rep*n.f)
      N.annual.jan.max  <- N.annual.jan.max/(n.rep*n.f)
      N.annual.jan.min  <- N.annual.jan.min/(n.rep*n.f)
      
      N.annual.dec      <- N.annual.dec/(n.rep*n.f)
      N.annual.dec.sd   <- N.annual.dec.sd/(n.rep*n.f)
      N.annual.dec.mean <- N.annual.dec.mean/(n.rep*n.f)
      N.annual.dec.max  <- N.annual.dec.max/(n.rep*n.f)
      N.annual.dec.min  <- N.annual.dec.min/(n.rep*n.f)
      
      # WRITE RESULTS
      
      if(writeTF == T){
        
        write(t(a.val),paste(dir,"/a_values",sep=""),ncolumns=n.a)
        write(t(b.val),paste(dir,"/b_values",sep=""),ncolumns=n.b)
        write(t(f.values),paste(dir,"/f_values",sep=""),ncolumns=n.f)
        write(t(N.annual.jan),paste(dir,"/n_annual_median",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.sd),paste(dir,"/n_annual_sd",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.mean),paste(dir,"/n_annual_mean",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.max),paste(dir,"/n_annual_max",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.min),paste(dir,"/n_annual_min",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.max.max),paste(dir,"/n_annual_max_max",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.min.min),paste(dir,"/n_annual_min_min",sep=""),ncolumns=n.f)
        
        write(t(N.annual.dec),paste(dir,"/n_annual_recr_median",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.sd),paste(dir,"/n_annual_recr_sd",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.mean),paste(dir,"/n_annual_recr_mean",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.max),paste(dir,"/n_annual_recr_max",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.min),paste(dir,"/n_annual_recr_min",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.max.max),paste(dir,"/n_annual_recr_max_max",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.min.min),paste(dir,"/n_annual_recr_min_min",sep=""),ncolumns=n.f)
        
        write(t(N.annual.jan.m),paste(dir,"/n_annual_m_median",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.sd),paste(dir,"/n_annual_m_sd",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.mean),paste(dir,"/n_annual_m_mean",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.max),paste(dir,"/n_annual_m_max",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.min),paste(dir,"/n_annual_m_min",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.max.max),paste(dir,"/n_annual_m_max_max",sep=""),ncolumns=n.f)
        write(t(N.annual.jan.m.min.min),paste(dir,"/n_annual_m_min_min",sep=""),ncolumns=n.f)
        
        write(t(N.annual.dec.m),paste(dir,"/n_annual_recr_m_median",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.sd),paste(dir,"/n_annual_recr_m_sd",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.mean),paste(dir,"/n_annual_recr_m_mean",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.max),paste(dir,"/n_annual_recr_m_max",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.min),paste(dir,"/n_annual_recr_m_min",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.max.max),paste(dir,"/n_annual_recr_m_max_max",sep=""),ncolumns=n.f)
        write(t(N.annual.dec.m.min.min),paste(dir,"/n_annual_recr_m_min_min",sep=""),ncolumns=n.f)
        
        write(t(N.killed),paste(dir,"/n_killed_median",sep=""),ncolumns=n.f)
        write(t(N.killed.sd),paste(dir,"/n_killed_sd",sep=""),ncolumns=n.f)
        write(t(N.killed.mean),paste(dir,"/n_killed_mean",sep=""),ncolumns=n.f)
        write(t(N.killed.max),paste(dir,"/n_killed_max",sep=""),ncolumns=n.f)
        write(t(N.killed.min),paste(dir,"/n_killed_min",sep=""),ncolumns=n.f)
        write(t(N.killed.max.max),paste(dir,"/n_killed_max_max",sep=""),ncolumns=n.f)
        write(t(N.killed.min.min),paste(dir,"/n_killed_min_min",sep=""),ncolumns=n.f)
        
        write(t(f.mort.12),paste(dir,"/f_mort_12",sep=""),ncolumns=n.f)
        write(t(prob.ext),paste(dir,"/prob_ext",sep=""),ncolumns=n.f)
        
        write(t(mean.t),paste(dir,"/mean_t",sep=""),ncolumns=n.f)
        write(t(median.t),paste(dir,"/median_t",sep=""),ncolumns=n.f)
        write(t(sd.t),paste(dir,"/sd_t",sep=""),ncolumns=n.f)
        write(t(max.t),paste(dir,"/max_t",sep=""),ncolumns=n.f)
        write(t(min.t),paste(dir,"/min_t",sep=""),ncolumns=n.f)
      }
      
  ###########################################################################################
  ##### PLOT SUMMARY PLOTS (4 plots) | auxiliary function 1 #################################
      
      plot.res <- function(f,dat,datmin,datmax,arrowTF,chosen,main){
        
        # select some graphics parameters
        
        colors <- c("red","blue","green","black")
        
        xlab <- expression(italic(f))
        logxy <- ""
        
        # plot a
        if(main == "a") {
          ylab="Number of animals 1st month recruitment"; logxy="xy"; xmin=min(f[f!=0]); xmax=1; 
          # ylab="Number of animals in May"; logxy="xy"; xmin=min(f[f!=0]); xmax=1; 
        }
        
        # plot b
        if(main == "b") {
          ylab="Number of animals killed per year"; logxy="x"; xmin=min(f[f!=0]); xmax=1; 
        }
        
        # plot d
        if(main == "d") {
          ylab="Time to extinction [years]"; pos="topleft"; 
          xmin <- 0; xmax <- 1
        }
        
        # draw plot
        
        ymax <- max(datmax[max(chosen),])
        res  <- datmin
        ymin <- min(res)
        
        #if(file=="n_annual_m_median") ymin <- 1 # plot a
        #if(file=="n_killed_median") ymin <- 0   # plot b
        
        res <- dat[1,]
        x <- f[res !=0]; y <- res[res !=0]
        
        plot(x, y, xlab=xlab, ylab=ylab, xlim=c(xmin,xmax), ylim=c(ymin,ymax),type="l", 
             main=main, log=logxy, col=colors[1], cex.main=1.5, cex.lab=1)
        
        for(i in 2:length(chosen)){
          
          res <- dat[chosen[i],]
          x <- f[res !=0]; y <- res[res !=0]
          
          lines(x, y, col=colors[i])
        }
        
        for(i in 1:length(chosen)){
          
          res <- dat[chosen[i],]
          x <- f[res !=0]; y <- res[res !=0]
          
          resmax <- datmax[chosen[i],]; resmax <- resmax[res != 0]
          resmin <- datmin[chosen[i],]; resmin <- resmin[res != 0]
          yl <- resmin
          yl[yl < 0] <- 0.001
          yt <- resmax
          
          # add arrows or not 
          
          if(arrowTF==T){
            arrows(x, yl, x, yt, length=0.05, angle=90, code=3, col="red")
            lines(x,y)
          } else {
            lines(x,yl,lty=3,col=colors[i])
            lines(x,yt,lty=2,col=colors[i])
          }
        }
        
        #if(file=="n_annual_m_median")
        legend("topright",c("max","median","min"),lty=c(2,1,3), bty="n", cex=0.65)
      }
      
  ###########################################################################################
  ##### PLOT SUMMARY PLOTS (4 plots) | auxiliary function 2 #################################
      
      plot.res.pe <- function(f,dat, median.pop.size, chosen){  # results probability of extinction
        
        # draw plot
        # select some graphics parameters
        
        main <- "c" 
        xlab <- expression(italic(f))
        ylab="Probability of extinction"
        pos <- "bottomright"
        cols <- c("red","blue","green","black")
        
        xmin <- 0; xmax <- 1
        ymax <- max(dat);
        ymin <- min(dat)
        if(ymax == ymin) {ymax <- 1; ymin <- 0}
        
        for(i in 1:length(chosen)){
          
          y <- dat[chosen[i],] + 0.0001
          x <- f; 
          
          if(i==1){
            
            plot(x, y, xlab=xlab, ylab=ylab, xlim=c(xmin,xmax), 
                 ylim=c(ymin,ymax),type="l", main=main, 
                 col=cols[1], cex.main=1.5, cex.lab=1)
          }
          else
            lines(x, y, lty=1, col=cols[i])
        }
        
        # now add legend
        
        txt <- paste("Median orig. pop. =",trunc(median.pop.size[chosen[1]]),sep="")
        if(length(chosen)>1)
          for(i in 2:length(chosen))
            txt <- c(txt,paste("Median orig. pop. =",trunc(median.pop.size[chosen[i]]),sep=""))
        
        legend(pos,txt,lty=1,col=cols, bty="n", cex=0.65)
      }
      
  ###########################################################################################
  ##### PLOT SUMMARY PLOTS (4 plots) | the plots ############################################
      
      if(plotTF ==T){
        
        arrowTF=F
        chosen.a <- c(4,5,6,7)
        
        # quartz()
        layout(matrix(1:4,ncol=2,byrow=T))
        
        a.val <- scan(paste(dir,"/a_values",sep=""))
        
        chosen <- which(a.val == chosen.a[1])
        if(length(chosen.a) > 1)
          for(i in 2:length(chosen.a))
            chosen <- c(chosen,which(a.val == chosen.a[i]))
        
        plot.res(f.values,N.annual.jan.m, N.annual.jan.m.min.min, N.annual.jan.m.max.max, 
                 arrowTF,chosen,"a")
        plot.res(f.values,N.killed, N.killed.min.min, N.killed.max.max, arrowTF=arrowTF,chosen,"b")
        plot.res.pe(f.values,prob.ext, N.annual.jan, chosen)
        plot.res(f.values,mean.t, min.t, max.t, arrowTF=arrowTF,chosen,"d")
        
      }
      
      cat("start.time = ", start.time, "\n")
      cat("end.time = ", date(), "\n", "\n")
      
  ###########################################################################################
  #-----------------------------------------------------------------------------------------#
  return(list("n.rep" = n.rep, "years" = years, 
              "transient1" = transient1, "transient2" = transient2, "transient3" = transient3, 
              "f.values" = f.values, "a.val" = a.val, "b.val" = b.val, # "b" = b, 
              "n.ages" = n.ages, "min.age" = min.age, "repr.period" = repr.period,
              "months.of.recruitment" = months.of.recruitment, 
              "prob.recr.month.m" = prob.recr.month.m, 
              "prob.recr.month.sd" = prob.recr.month.sd,  
              "prop.recr.age.m" = prop.recr.age.m, "prop.recr.age.sd" = prop.recr.age.sd, 
              "clutch.size.m" = clutch.size.m, "clutch.size.sd" = clutch.size.sd, 
              "jsmm" = jsmm, "jssdm" = jssdm, "jsmf" = jsmf, "jssdf" = jssdf, 
              "ysmm" = ysmm, "yssdm" = yssdm, "ysmf" = ysmf, "yssdf" = yssdf, 
              "asmm" = asmm, "assdm" = assdm, "asmf" = asmf, "assdf" = assdf,
              "plotTF" = plotTF, "writeTF" = writeTF, "dirname" = dirname.rep))    
  #-----------------------------------------------------------------------------------------#
  ###########################################################################################
  }
  # end of agemodel_rep
  #-------------------------------------------------------------------------------------#
  #######################################################################################
  #
  #
  #######################################################################################
  #-------------------------------------------------------------------------------------#
  # evaluate what input is given for the parameters
  if (isTRUE(identical(par, "")) && isTRUE(identical(file, "")) )
    {
      # if the user does not provide a list or a file with the parameter values
      # use default file
      file.i <- system.file("extdata", "pam.txt", package="PopIn")
      con.i <- file(file.i, "r")
      while ( TRUE ) {
        line <- readLines(con.i, n = 1)
        if ( length(line) == 0 ) {
          break
        }
        eval(parse(text = line))
      }
      close(con.i)
    }
    
    else if (isFALSE(identical(par, "")) && isFALSE(identical(file, "")) )
    { 
      # if the user provides a list AND a file with the parameter values
      # return warning message
      cat("too much information... use a list OR a file with the parameter values", "\n")
      return(NULL)
    }
    
    else if(any(param!="")){
      # if the user provides a list with the parameter values
      # read user list
      years <- param$years
      transient1 <- param$transient1
      transient2 <- param$transient2
      transient3 <- param$transient3
      f.mort <- param$f.mort
      a <- param$a
      b <- param$b
      f.values <- param$f.values
      a.val <- param$a.val
      n.ages <- param$n.ages
      min.age <- param$min.age
      repr.period <- param$repr.period
      months.of.recruitment <- param$months.of.recruitment
      prob.recr.month.m <- param$prob.recr.month.m
      prob.recr.month.sd <- param$prob.recr.month.sd
      prop.recr.age.m  <- param$prop.recr.age.m
      prop.recr.age.sd <- param$prop.recr.age.sd
      clutch.size.m <- param$clutch.size.m
      clutch.size.sd <- param$clutch.size.sd
      jsmm <- param$jsmm
      jssdm <- param$jssdm
      jsmf <- param$jsmf
      jssdf <- param$jssdf
      ysmm <- param$ysmm
      yssdm <- param$yssdm
      ysmf <- param$ysmf
      yssdf <- param$yssdf
      asmm <- param$asmm
      assdm <- param$assdm
      asmf <- param$asmf
      assdf <- param$assdf
      plotTF <- param$plotTF
      writeTF <- param$writeTF
      dirname.basic <- param$dirname.basic
      dirname.rep <- param$dirname.rep
    }
    
    else{
      # if the user provides a file with the parameter values
      # read user file
      con <- file(file, "r")
      while ( TRUE ) {
        line <- readLines(con, n = 1)
        if ( length(line) == 0 ) {
          break
        }
        eval(parse(text = line))
      }
      close(con)
    }
    
  #######################################################################################
  #-------------------------------------------------------------------------------------#
  # now evaluate the number of replicates 
  #  to choose which function will be called
    
    if(n.rep==1){
      # call agemodel_basic
      agemodel_basic(years,transient1, transient2, transient3, 
                     f.mort, a, b, 
                     n.ages, min.age, repr.period, months.of.recruitment, 
                     prob.recr.month.m, prob.recr.month.sd, 
                     prop.recr.age.m, prop.recr.age.sd, 
                     clutch.size.m, clutch.size.sd, 
                     jsmm, jssdm, jsmf, jssdf, 
                     ysmm, yssdm, ysmf, yssdf, 
                     asmm, assdm, asmf, assdf,
                     plotTF, writeTF, dirname.basic)
    }
    
    else{
      # if n.rep > 1
      # call agemodel_rep
      agemodel_rep(n.rep, years, transient1, transient2, transient3, 
                   a.val, f.values, 
                   n.ages, min.age, repr.period, months.of.recruitment, 
                   prob.recr.month.m, prob.recr.month.sd, 
                   prop.recr.age.m, prop.recr.age.sd, 
                   clutch.size.m, clutch.size.sd,
                   jsmm, jssdm, jsmf, jssdf, 
                   ysmm, yssdm, ysmf, yssdf, 
                   asmm, assdm, asmf, assdf,
                   plotTF, writeTF, dirname.rep)
    }
} 
##########################################################################################
#
#                    Recr              Recr
#                     |                 |
#              Z + RM | Z + RM   Z + RM | Z + RM
#             ------>\/-------> ------>\/------->
#    --------|--------|--------|--------|--------
#           May     June     July     August
#
#    N(June) = N(May)*exp(-Z) - Road.Mort(May) + Recr
#
#    Recr = f(N(May)*exp(-Z))
#
#    N(July) = N(June)*exp(-Z) - Road.Mort(June)
#
##################################### end of agemodel #####################################
#-----------------------------------------------------------------------------------------#
#
############# end of PopIn - impact of linear Infrastructures on Populations ##############
#-----------------------------------------------------------------------------------------#