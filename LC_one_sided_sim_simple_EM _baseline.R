setwd("C:\\UCLA\\thesis_ideas\\PhD_thesis\\latent_class_model\\LC_one_sided_sim_EM")
# source("LC_one_sided_sim_simple_EM.R")

rm(list = ls(all = TRUE))          

library("mlogit")
library("ggplot2")
# library(gmnl)

source("glogitform.R")
source("gmnl.draws.R")
source("gmnl.logliks.R")
source("gmnl.methods.R")
source("gmnl.tools.R")
source("gmnl.R")

set.seed(1234)               

# simulation settings
B = 3

# Globals
intercept = TRUE
outside_opt = TRUE
N_array <- c(400) # , 500, 600)      # Number of individuals
J <- 4                      # Number of alternatives
t_array <- c(30) # , 50 , 60)        # Number of choice situations
Q = 2                       # Number of classes

# b_truth = c(0.5, 0.5, 1.0,  # class 1: b1, b2, b3
#             1.0, 1.0, 0.5)  # class 2: b1, b2, b3
# num_beta = 3
# att_index = 4:6

b_truth = c(0.5, 0.5, 1.0, 1.0,  # class 1: b1, b2, b3, b4
            1.0, 1.0, 0.5, 0.5)  # class 2: b1, b2, b3, b4
num_beta = 4
att_index = 4:7

b_truth_mat=matrix(b_truth, nrow=num_beta)

# store simulation result
if (intercept) {
  if (outside_opt) {
    temp = (num_beta + J)*Q
  } else {
    temp = (num_beta + J-1)*Q
  }
} else {
  temp = num_beta*Q
}
params_all = matrix(0, ncol=temp+3+Q, nrow=B*length(t_array)*length(N_array))
params_all_idx = 1

temp = NULL
if (intercept) {
  if (outside_opt) {
    for (i in 1:Q) {
      temp <- c(temp, paste("class",i, 2:(J+1), "intercept",sep="."), paste("class", i, "beta", 1:num_beta, sep="."))
    }
  } else {
    for (i in 1:Q) {
      temp <- c(temp, paste("class",i, 2:J, "intercept",sep="."), paste("class", i, "beta", 1:num_beta, sep="."))
    }
  }
} else {
  for (i in 1:Q) {
    temp <- c(temp, paste("class", i, "beta", 1:num_beta, sep="."))
  }
}
colnames(params_all) = c(temp, "T", "N", "MaxLL", paste("class", 1:Q, "share", sep="."))

# has_NA = vector(mode="logical", length=B*length(t_array)*length(N_array))

for (N in N_array) {
  
  b1 = rep(b_truth_mat[1,], each=N/Q)
  b2 = rep(b_truth_mat[2,], each=N/Q)
  b3 = rep(b_truth_mat[3,], each=N/Q)
  b4 = rep(b_truth_mat[4,], each=N/Q)
  
  
  for (t in t_array) {
    
    for (b in 1:B) {
      
      print(paste0("N = ", N))
      print(paste0("t = ", t))
      print(paste0("b = ", b))
      
      id <- rep(1:N, each = J * t) # Id for each individual
      
      # Create random parameters and variables
      # x1 = rnorm(N * J * t)                   # continuous variable
      # x2 = as.numeric(runif(N * J * t) < 0.5) # dummy variable
      # # x2 = rnorm(N * J * t) 
      # x3 = rnorm(N * J * t)                   # continuous variable
      # x4 = rnorm(N * J * t)
      
      x1 = cut(runif(N * J * t), breaks = c(0, 0.2,0.5, 0.6, 1), labels = 1:4) 
      x1 = as.numeric(levels(x1))[x1]
      x2 = as.numeric(runif(N * J * t) < 0.5) # dummy variable
      x3 = rnorm(N * J * t)                   
      x3 = cut(x3, breaks = seq(min(x3), max(x3), length=4), labels=1:3, include.lowest=TRUE)
      x3 = as.numeric(levels(x3))[x3]
      x4 = as.numeric(runif(N * J * t) < 0.6) # dummy variable
      
      
      theta_0 = rep(0,num_beta*Q)
      # theta_0 = c(0.7313592, 0.7566393, 0.6730434,
      #             0.7713592, 0.7966393, 0.7130434)
      # theta_0 = c(0.7097857, 0.7317511, 0.6411065,
      #             0.7497857, 0.7717511, 0.6811065)
      
      # The true data generating process
      eta <- -log(-log(runif(N * J * t)))
      
      # adjust for outside option (not buying anything)
      if (outside_opt) {
        eta0 = rep(apply(-log(-log(matrix(runif(N * J * t), J))), 2, max), each=J)
        eta = eta - eta0
      }
      
      # U <-  b1[id] * x1 + b2[id] * x2 + b3[id] * x3 + eta
      U <-  b1[id] * x1 + b2[id] * x2 + b3[id] * x3 + b4[id] * x4 + eta
      
      # Chosen alternative
      choice <- rep(0, N * J * t)
      for (i in 1:(N * t)) {
        U_j <- U[((i - 1) * J + 1):(i * J)]
        U_max <- max(U_j)
        if (!outside_opt || (U_max > 0)) {
          pos <- which(U_j == U_max)
          choice[((i - 1) * J + 1):(i * J)][pos] <- TRUE
        }
      }
      
      data <- as.data.frame(cbind(id, choice, alt=rep(1:J, N*t), x1, x2, x3, x4))
      # data <- as.data.frame(cbind(id, choice, alt=rep(1:J, N*t), x1, x2, x3))
      data <- mlogit.data(data, 
                          choice = "choice",
                          shape = "long", 
                          alt.levels = c("1", "2", "3", "4"), 
                          id.var = "id")
      
      ##################################### latent class model ####################################
      # pad the data with the outside option (so we can use the gmnl function withoud modification)
      if (outside_opt) {
        temp = lapply(split(data, rep(1:(N*t), each=J)), function(x) rbind(x,0))
        data_new = Reduce(rbind, temp)
        rownames(data_new)[seq(from=(J+1),to=nrow(data_new),by=(J+1))]= paste(1:(N*t), (J+1), sep=".")
        data_new$id[seq(from=(J+1),to=nrow(data_new),by=(J+1))] = data_new$id[(seq(from=(J+1),to=nrow(data_new),by=(J+1)))-1]
        data_new$alt[seq(from=(J+1),to=nrow(data_new),by=(J+1))] = (J+1)
        new_choice = Reduce(rbind, split(choice, rep(1:J, times=N*t/J)))
        temp = apply(new_choice, 2, function(x) all(x==0))
        data_new$choice[seq(from=(J+1),to=nrow(data_new),by=(J+1))] = as.logical(temp)
        
        data_new <- mlogit.data(data_new, 
                                choice = "choice",
                                shape = "long", 
                                alt.levels = c("1", "2", "3", "4", "5"), 
                                id.var = "id")
        data = data_new
      }
      
      if (intercept) {
        out.lc <- gmnl(choice ~ x1 + x2 + x3 + x4 | 1 | 0 | 0 | 1, print.init=TRUE, 
                           data = data, model = 'lc', panel = TRUE, Q = Q)
      } else {
        out.lc <- gmnl(choice ~ x1 + x2 + x3 + x4 | 0 | 0 | 0 | 1, print.init=TRUE,
                       data = data, model = 'lc', panel = TRUE, Q = Q)
      }
      
      
      
      # out.lc <- gmnl_debug(choice ~ x1 + x2 + x3 + x4 | 0 | 0 | 0 | 1, print.init=TRUE, start = theta_0,
      #                      data = data, model = 'lc', panel = TRUE, Q = Q, att_index=att_index, outside_opt=outside_opt,
      #                      opt_algo="bfgs")
      #
      
      ############################################################################################
      
      print(summary(out.lc))
      
      # get the coefficients
      if (intercept) {
        if (outside_opt) {
          coef_lc = coef(out.lc)[c((J+1):(J+num_beta), (J+num_beta+J+1):(J+num_beta+J+num_beta))]
        } else {
          coef_lc = coef(out.lc)[c(J:(J+num_beta-1), (J+num_beta+J-1):(J+num_beta+J-1+num_beta-1))]
        }
      } else {
        coef_lc = coef(out.lc)
      }
      
      # shares of classes
      class_idx = length(coef(out.lc))-(Q-2)
      class_total = exp(0) + sum(exp(coef(out.lc)[class_idx:length(coef(out.lc))]))
      class_shares = c(exp(0), exp(coef(out.lc)[class_idx:length(coef(out.lc))]))/class_total
      print(paste0("Class ", 1:Q, " share: ", class_shares))
      
      # the classes may not be in the same order as expected
      # so i'm trying to rearrange them here (currently only works for Q=2)
      dist1 = sum((b_truth - coef_lc[1:(num_beta*Q)])^2)
      dist2 = sum((b_truth - coef_lc[c((num_beta+1):(num_beta*2), 1:num_beta)])^2)
      
      if (is.null(out.lc$maximum)) out.lc$maximum = NA
      
      if (dist1 <= dist2)
      {
        params_all[params_all_idx,] = c(coef(out.lc)[1:(class_idx-1)], t, N, out.lc$maximum, class_shares)
      } else {
        if (intercept) {
          temp = coef(out.lc)[1:(class_idx-1)]
          params_all[params_all_idx,] = c(coef(out.lc)[c((length(temp)/Q + 1):length(temp), 1:(length(temp)/Q))],
                                          t, N, out.lc$maximum, class_shares[Q:1])
        } else {
          params_all[params_all_idx,] = c(coef_lc[c((num_beta+1):(num_beta*2), 1:num_beta)], 
                                          t, N, out.lc$maximum, class_shares[Q:1])
        }
        # just to be consistent
        coef_lc = coef_lc[c((num_beta+1):(num_beta*2), 1:num_beta)]
      }
      params_all_idx = params_all_idx + 1
      
    } # end of b loop
    
  } # end of t loop
} # end of N loop

# plot truth
hline.dat = data.frame(beta=paste0("beta.",rep(1:num_beta,times=Q)), coeff=b_truth, class=rep(paste("class",1:Q,sep="."), each=num_beta))

if (intercept) {
  if (outside_opt) {
    if (outside_opt) {
      coef_lc_all = params_all[, c((J+1):(J+num_beta), (J+num_beta+J+1):(J+num_beta+J+num_beta))]
    } else {
      coef_lc_all = params_all[, c(J:(J+num_beta-1), (J+num_beta+J-1):(J+num_beta+J-1+num_beta-1))]
    }
  }
  
  params_all_df = data.frame(coeff = matrix(coef_lc_all, ncol=1),
                             beta = as.factor(rep(rep(paste("beta",1:num_beta,sep="."), each=B*length(t_array)*length(N_array)), times=Q)),
                             class = as.factor(rep(paste("class", 1:Q, sep="."), each=num_beta*B*length(t_array)*length(N_array))),
                             t = as.factor(rep(rep(t_array, each=B, times=length(N_array)), times=num_beta*Q)),
                             N = as.factor(rep(rep(N_array, each=B*length(t_array)), times=num_beta*Q)))
} else {
  params_all_df = data.frame(coeff = matrix(params_all[,1:(num_beta*Q)], ncol=1),
                             beta = as.factor(rep(rep(paste("beta",1:num_beta,sep="."), each=B*length(t_array)*length(N_array)), times=Q)),
                             class = as.factor(rep(paste("class", 1:Q, sep="."), each=num_beta*B*length(t_array)*length(N_array))),
                             t = as.factor(rep(rep(t_array, each=B, times=length(N_array)), times=num_beta*Q)),
                             N = as.factor(rep(rep(N_array, each=B*length(t_array)), times=num_beta*Q)))
}
 
  
  


for (nn in N_array) {
  
  print(ggplot(data=subset(params_all_df, N==nn), aes(x=t, y=coeff)) + geom_boxplot() + facet_grid(class ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =coeff), color="red") + 
          ggtitle(paste0("intercept = ", intercept, ", outside option = ", outside_opt, ", N = ", nn, ", B = ", B)))
  
  png(paste0("LC_int_", intercept, "_outside_opt_", outside_opt, "_B" ,B, "_N", nn, ".png"))
  print(ggplot(data=subset(params_all_df, N==nn), aes(x=t, y=coeff)) + geom_boxplot() + facet_grid(class ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =coeff), color="red") + 
          ggtitle(paste0("intercept = ", intercept, ", outside option = ", outside_opt, ", N = ", nn, ", B = ", B)))
  dev.off()
}

for (nn in N_array) {
  for (tt in t_array) {
    df = as.data.frame(params_all[which(params_all[,"N"]==nn & params_all[,"T"]==tt),])
    beta_med = apply(df[,1:(num_beta*Q)],2,median)
    print(paste0("N = ", nn, " T = ", tt, ": "))
    print(beta_med)
  }
}
