#' Power and Sample Size Estimation for the Win Ratio
#' @description A simulation-based power and sample size calculation approach, relying on the relationship between the win ratio and the rank distribution
#'
#' @param n_arm_1 The sample size in arm 1
#' @param n_arm_2 The sample size in arm 2
#' @param alpha The two-sided type I error rate
#' @param WinRatio1 The specified Win Ratio for ratio of Arm 1 wins to Arm 1 losses on the first endpoint, either p1 or WinRatio1 should be specified
#' @param WinRatio2 The specified Win Ratio for ratio of Arm 1 wins to Arm 1 losses on the second endpoint, either p2 or WinRatio2 can be specified. Do not include if only assuming one endpoint
#' @param p1 The probability of success in arm 1 for the first endpoint, either p1 or WinRatio1 should be specified
#' @param p2 The probability of success in arm 1 for the second endpoint, either p2 or WinRatio2 can be specified. Do not include if only assuming one endpoint
#' @param corr The correlation of endpoints. Must be specified if p2 or WinRatio2 is specified.
#' @param cens_prop The proportion of administrative censoring on the first endpoint
#' @param ties The proportion of ties on the second endpoint
#' @param ties_dist Equal to either "best" if the better ranks are tied or "worst" if worse ranks are tied
#' @param n.iter The number of simulations
#' @param seed The seed for random number generation
#' @param boot The number of bootstrap samples within each simulation
#' @param print.iter If TRUE, each simulation number will be printed after completion
#' @return A list with components including the median, 25th percentile, 75th percentile, mean, and standard deviation of the true win ratios assuming complete data and the
#' observed win ratios assuming censoring and/or ties across simulations;
#' and the power of the proposed sample size under specified assumptions
#' @examples sim_example<-WinRatio_power(n_arm_1=10,n_arm_2=10,alpha=0.05,p1=0.6,n.iter=100,seed=1234,boot=100)
#' @export

WinRatio_power <- function(n_arm_1, n_arm_2, alpha=0.05, WinRatio1=NULL, p1=NULL, WinRatio2=NULL, p2=NULL, corr=0,
                              cens_prop=0, ties=0, ties_dist="best", n.iter=500, seed=1000, boot=500, print.iter=FALSE){

  #Warnings
  if (!is.null(alpha) && !is.numeric(alpha) || any(alpha < 0 | alpha > 1))
    stop("'alpha' must be numeric and in [0, 1]")

  if (!is.null(n_arm_1) && !is.numeric(n_arm_1) || any(n_arm_1 <= 0))
    stop("'n_arm_1' must be numeric and > 0")

  if (!is.null(n_arm_2) && !is.numeric(n_arm_2) || any(n_arm_2 <= 0))
    stop("'n_arm_2' must be numeric and > 0")

  if (!is.null(WinRatio1) && !is.numeric(WinRatio1) || any(WinRatio1 <= 0))
    stop("'WinRatio1' must be numeric and > 0")

  if (!is.null(p1) && !is.numeric(p1) || any(p1 <= 0 | p1 >= 1))
    stop("'p1' must be numeric and in [0, 1]")

  if (!is.null(WinRatio2) && !is.numeric(WinRatio2) || any(WinRatio2 <= 0))
    stop("'WinRatio2' must be numeric and > 0")

  if (!is.null(p2) && !is.numeric(p2) || any(p2 <= 0 | p2 >= 1))
    stop("'p2' must be numeric and in [0, 1]")

  if (is.null(p1) & is.null(WinRatio1))
    stop("either 'p1' or 'WinRatio1' needs to be specified")

  if(!is.null(p1) & !is.null(WinRatio1) & isTRUE(p1!=(WinRatio1/(WinRatio1+1))))
    stop("'p1' and 'WinRatio1' are both specified but p1 not equal to WinRatio1/(WinRatio1+1)")

  if(!is.null(ties) & is.null(ties_dist))
    stop("`ties_dist` must be specified if `ties` > 0")

  #If WinRatio1 not specified, use "p1"
  if (is.null(WinRatio1)){p1<-p1}
  #If WinRatio1 specified, convert to "p1"
  if (is.null(p1)){p1<-WinRatio1/(WinRatio1+1)}

  #If WinRatio2 not specified, but p2 is specified
  if (is.null(WinRatio2) & !is.null(p2)){p2<-p2}
  #If WinRatio2 specified, convert to "p2"
  if (is.null(p2)){p2<-WinRatio2/(WinRatio2+1)}

  no_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(no_cores)
  registerDoParallel(cl)

  n_total<-n_arm_1+n_arm_2
  odd<-ifelse(n_total%%2 == 1,1,0)


  wr_true<-rep(NA,n.iter)
  wr_est<-rep(NA,n.iter)
  wr_2<-rep(NA,n.iter)
  p_tie_observed<-rep(NA,n.iter)
  lower_ci_bca<-rep(NA,n.iter)
  upper_ci_bca<-rep(NA,n.iter)
  power_bca<-rep(NA,n.iter)

  for (j in 1:n.iter){
    set.seed(j+seed)

## If no censoring, only need first endpoint:
  if (cens_prop == 0){

    if (p1==0.5){
      #hypergeometric - under null
      arm_a_n_top_p1<-rhyper(1,0.5*(n_arm_1+n_arm_2),0.5*(n_arm_1+n_arm_2),m_arm_1)
    }

    # NEED TO CONFIRM THIS IS ACCURATE FOR RATIO != 1 (????)

    if (p1!=0.5){
      #noncentral hypergeometric - under alternative
      if (n_arm_1/n_arm_2 == 1){
      arm_1_n_top_p1 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p1)/log(p1))
      }
      if (n_arm_1/n_arm_2 != 1){
      x<- -n_arm_1*0.5 + n_arm_1**2/(n_arm_1+n_arm_2) + (p1/(n_arm_1+n_arm_2))*2*n_arm_1*((n_arm_1+n_arm_2)-n_arm_1)
      p1_avg<-x/n_arm_1
      arm_1_n_top_p1 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p1_avg)/log(p1_avg))
      }
    }

    if (odd==0){
      n_upper_rank<-n_total*0.5
      n_upper_rank_1<-n_upper_rank+1
    }
    if (odd==1){
      set.seed(seed+j+1)
      n_upper_rank<-floor(n_total*0.5) + rbinom(1,1,0.5)
      n_upper_rank_1<-n_upper_rank+1
    }

    arm_1_ranks_p1 <- c(sample(1:n_upper_rank,arm_1_n_top_p1), sample(n_upper_rank_1:(n_arm_1+n_arm_2),n_arm_1-arm_1_n_top_p1,))
    arm_1_ranks_scramble_p1 <- sample(arm_1_ranks_p1,n_arm_1)


    #assign remaining ranks to arm b
    arm_2_ranks_p1 <- setdiff(1:(n_arm_1+n_arm_2),arm_1_ranks_p1)
    arm_2_ranks_scramble_p1 <- sample(arm_2_ranks_p1,n_arm_2)

    #make adjacency matrix
    adj_mat <- matrix(0,n_arm_1+n_arm_2,n_arm_1+n_arm_2)
    adj_mat <- 1*upper.tri(adj_mat)
    adj_mat_o <- adj_mat #save original adjacency matrix

    arm_1_n_wins_true <- sum(rowSums(adj_mat_o[arm_1_ranks_p1,arm_2_ranks_p1]))
    arm_2_n_wins_true <- sum(rowSums(adj_mat_o[arm_2_ranks_p1,arm_1_ranks_p1]))
    wr_true[j]<-arm_1_n_wins_true/arm_2_n_wins_true

    wr_est[j] <- wr_true[j]





    #Bootstrap resample
    #estimate of win ratio for each bootstrap sample within one iteration
    winratio_est<-foreach::foreach(i=1:boot, .combine=c) %dopar% {
      set.seed(10000*j + i + seed)
      sample_treat<-sample(arm_1_ranks_p1,replace=T,size=n_arm_1)
      sample_control<-sample(arm_2_ranks_p1,replace=T,size=n_arm_2)
      arm_1_n_wins_boot <- sum(rowSums(adj_mat_o[sample_treat,sample_control]))
      arm_2_n_wins_boot <- sum(rowSums(adj_mat_o[sample_control,sample_treat]))
      arm_1_n_wins_boot/arm_2_n_wins_boot
    }

    #BCA confidence interval
    #https://www.frontiersin.org/articles/10.3389/fpsyg.2019.02215/full
    #proportion of bootstrap estimates less than the original estimate
    proportion<-length(winratio_est[winratio_est<wr_est[j]])/boot
    z_0<-qnorm(proportion)
    #jackknife estimate for acceleration parameter
    ranks<-c(arm_1_ranks_p1,arm_2_ranks_p1)
    arm<-c(rep(1,n_arm_1),rep(0,n_arm_2))
    ranks_jackknife<-rep(NA,length(ranks))
    arm_jackknife<-rep(NA,length(ranks))
    wr_est_jk<-rep(NA,length(ranks))
    for (k in 1:length(ranks)){
      ranks_jackknife<-ranks[-k]
      arm_jackknife<-arm[-k]

      arm_1_ranks_p1_jk <- ranks_jackknife[arm_jackknife==1]
      arm_2_ranks_p1_jk <- ranks_jackknife[arm_jackknife==0]
      arm_1_n_wins_jk <- sum(rowSums(adj_mat_o[arm_1_ranks_p1_jk,arm_2_ranks_p1_jk]))
      arm_2_n_wins_jk <- sum(rowSums(adj_mat_o[arm_2_ranks_p1_jk,arm_1_ranks_p1_jk]))

      #calculate the observed WR estimate - jackknife estimate
      wr_est_jk[k] <- arm_1_n_wins_jk/arm_2_n_wins_jk
    }
    mean_jk<-mean(wr_est_jk)
    jk<-data.frame(cbind(wr_est_jk,mean_jk))
    jk$num<-(jk$mean_jk-jk$wr_est_jk)**3
    jk$den<-(jk$mean_jk-jk$wr_est_jk)**2
    a<-sum(jk$num)/(6*(sum(jk$den)**(3/2)))
    alpha_1<- pnorm(z_0 + ((z_0+qnorm(alpha/2))/(1-a*(z_0+qnorm(alpha/2)))))
    alpha_2<- pnorm(z_0 + ((z_0+qnorm(1-(alpha/2)))/(1-a*(z_0+qnorm(1-(alpha/2))))))
    winratio_est_order<-winratio_est[order(winratio_est)]

    lower_ci_bca<-quantile(winratio_est_order,alpha_1,na.rm=TRUE)
    upper_ci_bca<-quantile(winratio_est_order,alpha_2,na.rm=TRUE)

    power_bca<-ifelse(upper_ci_bca < 1 | lower_ci_bca > 1, 1, 0)

    if(print.iter==TRUE){print(j)}
  }

  ## If censoring > 0, need both endpoints:
  if (cens_prop > 0){
    if (p1==0.5){
      #hypergeometric - under null
      arm_1_n_top_p1<-BiasedUrn::rhyper(1,0.5*(n_arm_1+n_arm_2),0.5*(n_arm_1+n_arm_2),n_arm_1)
    }


    if (p1!=0.5){
      #noncentral hypergeometric - under alternative
      #select ranks for arm a - first endpoint
      if (n_arm_1/n_arm_2 == 1){
        arm_1_n_top_p1 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p1)/log(p1))
      }
      if (n_arm_1/n_arm_2 != 1){
        x<- -n_arm_1*0.5 + n_arm_1**2/(n_arm_1+n_arm_2) + (p1/(n_arm_1+n_arm_2))*2*n_arm_1*((n_arm_1+n_arm_2)-n_arm_1)
        p1_avg<-x/arm_a_n
        arm_1_n_top_p1 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p1_avg)/log(p1_avg))
      }
    }

    if (odd==0){
      n_upper_rank<-n_total*0.5
      n_upper_rank_1<-n_upper_rank+1
    }
    if (odd==1){
      set.seed(seed+j+1)
      n_upper_rank<-floor(n_total*0.5) + rbinom(1,1,0.5)
      n_upper_rank_1<-n_upper_rank+1
    }

    arm_1_ranks_p1 <- c(sample(1:n_upper_rank,arm_1_n_top_p1), sample(n_upper_rank_1:(n_arm_1+n_arm_2),n_arm_1-arm_1_n_top_p1,))
    arm_1_ranks_scramble_p1 <- sample(arm_1_ranks_p1,n_arm_1)

    #assign remaining ranks to arm b
    arm_2_ranks_p1 <- setdiff(1:(n_arm_1+n_arm_2),arm_1_ranks_p1)
    arm_2_ranks_scramble_p1 <- sample(arm_2_ranks_p1,n_arm_2)

    #make adjacency matrix
    adj_mat <- matrix(0,n_arm_1+n_arm_2,n_arm_1+n_arm_2)
    adj_mat <- 1*upper.tri(adj_mat)
    adj_mat_o <- adj_mat #save original adjacency matrix

    arm_1_n_wins_true <- sum(rowSums(adj_mat_o[arm_1_ranks_p1,arm_2_ranks_p1]))
    arm_2_n_wins_true <- sum(rowSums(adj_mat_o[arm_2_ranks_p1,arm_1_ranks_p1]))
    wr_true[j]<-arm_1_n_wins_true/arm_2_n_wins_true

    ##effect size on second endpoint
    if (p2==0.5){
      #hypergeometric - under null
      arm_1_n_top_p2<-BiasedUrn::rhyper(1,0.5*(n_arm_1+n_arm_2),0.5*(n_arm_1+n_arm_2),n_arm_1)
    }
    if (p2!=0.5){
      #non-central hypergeometric - under alternative
      if (n_arm_1/n_arm_2 == 1){
        arm_1_n_top_p2 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p2)/log(p2))
      }
      if (n_arm_1/n_arm_2 != 1){
        x<- -n_arm_1*0.5 + n_arm_1**2/(n_arm_1+n_arm_2) + (p2/(n_arm_1+n_arm_2))*2*n_arm_1*((n_arm_1+n_arm_2)-n_arm_1)
        p2_avg<-x/n_arm_1
        arm_1_n_top_p2 <- BiasedUrn::rWNCHypergeo(nran=1,m1=n_arm_1,m2=n_arm_2,n=n_arm_1, odds=log(1-p2_avg)/log(p2_avg))
      }
    }
    if (odd==0){
      n_upper_rank<-n_total*0.5
      n_upper_rank_1<-n_upper_rank+1
    }
    if (odd==1){
      set.seed(seed+j+1)
      n_upper_rank<-floor(n_total*0.5) + rbinom(1,1,0.5)
      n_upper_rank_1<-n_upper_rank+1
    }

    arm_1_ranks_p2 <- c(sample(1:n_upper_rank,arm_1_n_top_p2), sample(n_upper_rank_1:(n_arm_1+n_arm_2),n_arm_1-arm_1_n_top_p2,))
    arm_1_ranks_scramble_p2 <- sample(arm_1_ranks_p2,n_arm_1)

    arm_2_ranks_p2 <- setdiff(1:(n_arm_1+n_arm_2),arm_1_ranks_p2)
    arm_2_ranks_scramble_p2 <- sample(arm_2_ranks_p2,n_arm_2)


    #wr - second endpoint
    arm_1_n_wins_p2 <- sum(rowSums(adj_mat_o[arm_1_ranks_p2,arm_2_ranks_p2]))
    arm_2_n_wins_p2 <- sum(rowSums(adj_mat_o[arm_2_ranks_p2,arm_1_ranks_p2]))
    wr_2[j]<-arm_1_n_wins_p2/arm_2_n_wins_p2


    #try the copula trick to induce dependence
    set.seed(10000+j)
    rc <- copula::rCopula(n=n_arm_1,copula=copula::normalCopula(corr))

    arm_1_dep_p1 <- rep(NA,n_arm_1)
    arm_1_dep_p2 <- rep(NA,n_arm_1)

    arm_1_dep_p1[order(rc[,1])] <- sort(arm_1_ranks_p1)
    arm_1_dep_p2[order(rc[,2])] <- sort(arm_1_ranks_p2)

    correlation_1<-cor(arm_1_dep_p1,arm_1_dep_p2,method="spearman")

    #Arm B
    set.seed(20000+j)
    rc_2 <- copula::rCopula(n=n_arm_2,copula=copula::normalCopula(corr))

    arm_2_dep_p1 <- rep(NA,n_arm_2)
    arm_2_dep_p2 <- rep(NA,n_arm_2)

    arm_2_dep_p1[order(rc_2[,1])] <- sort(arm_2_ranks_p1)
    arm_2_dep_p2[order(rc_2[,2])] <- sort(arm_2_ranks_p2)

    correlation_2<-cor(arm_2_dep_p1,arm_2_dep_p2,method="spearman")

    ##If ties on second endpoint
    if (ties > 0){

    tied_num<-rbinom(1,(n_arm_1+n_arm_2),ties)

    if (ties_dist == "worst"){
    ## Ties - worst ranks tied
    arm_1_dep_p2 <- ifelse(arm_1_dep_p2 >= (n_arm_1+n_arm_2)-tied_num, (n_arm_1+n_arm_2)-tied_num, arm_1_dep_p2)
    arm_2_dep_p2 <- ifelse(arm_2_dep_p2 >= (n_arm_1+n_arm_2)-tied_num, (n_arm_1+n_arm_2)-tied_num, arm_2_dep_p2)
    }

    if (ties_dist == "best"){
    ## Ties - best ranks tied
    arm_1_dep_p2 <- ifelse(arm_1_dep_p2 <= tied_num, 1, arm_1_dep_p2)
    arm_2_dep_p2 <- ifelse(arm_2_dep_p2 <= tied_num, 1, arm_2_dep_p2)
    }


    #wr - second endpoint removing ties
    arm_1_n_wins_p2_ties <- sum(rowSums(adj_mat_o[arm_1_dep_p2,arm_2_dep_p2]))
    arm_2_n_wins_p2_ties <- sum(rowSums(adj_mat_o[arm_2_dep_p2,arm_1_dep_p2]))
    wr_2_ties<-arm_1_n_wins_p2_ties/arm_2_n_wins_p2_ties

    p_tie_end2<- ((n_arm_1*n_arm_2)-(arm_1_n_wins_p2_ties+arm_2_n_wins_p2_ties))/(n_arm_1*n_arm_2)
    }


    data<-merge(arm_1_dep_p1,arm_2_dep_p1,all=TRUE)
    data_p2<-merge(arm_1_dep_p2,arm_2_dep_p2,all=TRUE)
    colnames(data_p2)[1:2]<-c("x_p2","y_p2")
    data<-cbind(data,data_p2)

    ##Assuming administrative censoring
    censored<-rbinom(1,(n_arm_1+n_arm_2),cens_prop)
    arm_1_censored<-ifelse(arm_1_dep_p1<= censored,1,0)
    arm_2_censored<-ifelse(arm_2_dep_p1<= censored,1,0)

    censoring <- (length(arm_1_censored[arm_1_censored==1]) +length(arm_2_censored[arm_2_censored==1])) / (n_arm_1+n_arm_2)
    censoring_1<- length(arm_1_censored[arm_1_censored==1]) / n_arm_1
    censoring_2<- length(arm_2_censored[arm_2_censored==1]) / n_arm_2

    data_1<-data.frame(cbind(arm_1_dep_p1,arm_1_censored))
    colnames(data_1)[1]<-c("x")
    data_2<-data.frame(cbind(arm_2_dep_p1,arm_2_censored))
    colnames(data_2)[1]<-c("y")
    data<-merge(data,data_2,by=c("y"),all.x=TRUE)
    data<-merge(data,data_1,by=c("x"),all.x=TRUE)

    #flip-able edges
    data$flippable<-ifelse(data$arm_1_censored==1 & data$arm_2_censored==1, 1, 0)

    if (length(data$flippable[data$flippable==1]) > 0){

      flippable<-data[data$flippable==1,]
      flippable$true_win<-ifelse(flippable$x < flippable$y,1,0)

      #wins flip to losses
      flippable$win_flip<-ifelse(flippable$true_win==1 & (flippable$y_p2 < flippable$x_p2),1,0)
      edges_flipped_win<-flippable[flippable$win_flip==1,]

      #losses flip to wins
      flippable$loss_flip<-ifelse(flippable$true_win==0 & (flippable$y_p2 > flippable$x_p2),1,0)
      edges_flipped_loss<-flippable[flippable$loss_flip==1,]

      #end in tie
      flippable$tie<-ifelse(flippable$y_p2 == flippable$x_p2,1,0)
      edges_tied<-flippable[flippable$tie==1,]

      #Flip edges
      for (k in 1:nrow(edges_flipped_loss)){
        adj_mat[unlist(edges_flipped_loss[k,1]), unlist(edges_flipped_loss[k,2])]<-1
        adj_mat[unlist(edges_flipped_loss[k,2]), unlist(edges_flipped_loss[k,1])]<-0
      }
      for (k in 1:nrow(edges_flipped_win)){
        adj_mat[unlist(edges_flipped_win[k,1]), unlist(edges_flipped_win[k,2])]<-0
        adj_mat[unlist(edges_flipped_win[k,2]), unlist(edges_flipped_win[k,1])]<-1
      }
      #Edges -> Ties
      for (k in 1:nrow(edges_tied)){
        adj_mat[unlist(edges_tied[k,1]), unlist(edges_tied[k,2])]<-0
        adj_mat[unlist(edges_tied[k,2]), unlist(edges_tied[k,1])]<-0
      }
    }

    #use adjacency matrices to calculate the observed number of wins in each arm
    arm_1_n_wins <- sum(rowSums(adj_mat[arm_1_dep_p1,arm_2_dep_p1]))
    arm_2_n_wins <- sum(rowSums(adj_mat[arm_2_dep_p1,arm_1_dep_p1]))

    #calculate the observed WR estimate
    wr_est[j] <- arm_1_n_wins/arm_2_n_wins

    p_tie_observed[j] <- 1-((arm_1_n_wins+arm_2_n_wins)/(n_arm_1*n_arm_2))

    #Bootstrap resample
    #estimate of win ratio for each bootstrap sample within one iteration
    winratio_est<-foreach::foreach(i=1:boot, .combine=c) %dopar% {
      set.seed(10000*j + i + seed)
      sample_treat<-sample(arm_1_dep_p1,replace=T,size=n_arm_1)
      sample_control<-sample(arm_2_dep_p1,replace=T,size=n_arm_2)
      arm_1_n_wins_boot <- sum(rowSums(adj_mat[sample_treat,sample_control]))
      arm_2_n_wins_boot <- sum(rowSums(adj_mat[sample_control,sample_treat]))
      arm_1_n_wins_boot/arm_2_n_wins_boot
    }

    #BCA confidence interval
    #https://www.frontiersin.org/articles/10.3389/fpsyg.2019.02215/full
    #proportion of bootstrap estimates less than the original estimate
    proportion<-length(winratio_est[winratio_est<wr_est[j]])/boot
    z_0<-qnorm(proportion)
    #jackknife estimate for acceleration parameter
    ranks<-c(arm_1_dep_p1,arm_2_dep_p1)
    arm<-c(rep(1,n_arm_1),rep(0,n_arm_2))
    ranks_jackknife<-rep(NA,length(ranks))
    arm_jackknife<-rep(NA,length(ranks))
    wr_est_jk<-rep(NA,length(ranks))
    for (k in 1:length(ranks)){
      ranks_jackknife<-ranks[-k]
      arm_jackknife<-arm[-k]

      arm_1_dep_p1_jk <- ranks_jackknife[arm_jackknife==1]
      arm_2_dep_p1_jk <- ranks_jackknife[arm_jackknife==0]
      arm_1_n_wins_jk <- sum(rowSums(adj_mat[arm_1_dep_p1_jk,arm_2_dep_p1_jk]))
      arm_2_n_wins_jk <- sum(rowSums(adj_mat[arm_2_dep_p1_jk,arm_1_dep_p1_jk]))

      #calculate the observed WR estimate - jackknife estimate
      wr_est_jk[k] <- arm_1_n_wins_jk/arm_2_n_wins_jk
    }
    mean_jk<-mean(wr_est_jk)
    jk<-data.frame(cbind(wr_est_jk,mean_jk))
    jk$num<-(jk$mean_jk-jk$wr_est_jk)**3
    jk$den<-(jk$mean_jk-jk$wr_est_jk)**2
    a<-sum(jk$num)/(6*(sum(jk$den)**(3/2)))
    alpha_1<- pnorm(z_0 + ((z_0+qnorm(alpha/2))/(1-a*(z_0+qnorm(alpha/2)))))
    alpha_2<- pnorm(z_0 + ((z_0+qnorm(1-(alpha/2)))/(1-a*(z_0+qnorm(1-(alpha/2))))))
    winratio_est_order<-winratio_est[order(winratio_est)]

    lower_ci_bca[j]<-quantile(winratio_est_order,alpha_1,na.rm=TRUE)
    upper_ci_bca[j]<-quantile(winratio_est_order,alpha_2,na.rm=TRUE)

    power_bca[j]<-ifelse(upper_ci_bca[j] < 1 | lower_ci_bca[j] > 1, 1, 0)

    if(print.iter==TRUE){print(j)}
  }

}
    closeAllConnections()

  if(quantile(wr_est,1.0, na.rm=TRUE)=="Inf")
    warning("Some iterations producing infinity estimates for win ratio, consider larger sample size or less extreme win ratio")

  Power<-length(power_bca[power_bca==1])/length(power_bca)
  CI_width<-upper_ci_bca - lower_ci_bca

  mylist<-list("WR_True_Median" = summary(wr_true)[3], "WR_True_25th" = summary(wr_true)[2], "WR_True_75th"=summary(wr_true)[5],
               "WR_True_Mean" = summary(wr_true)[4], "WR_True_SD" = sd(wr_true),

               "WR_True2_Median" = summary(wr_2)[3], "WR_True2_25th" = summary(wr_2)[2], "WR_True2_75th"=summary(wr_2)[5],
               "WR_True2_Mean" = summary(wr_2)[4], "WR_True2_SD" = sd(wr_2),

               "WR_Obs_Median" = summary(wr_est)[3], "WR_Obs_25th" = summary(wr_est)[2], "WR_Obs_75th"=summary(wr_est)[5],
               "WR_Obs_Mean" = summary(wr_est)[4], "WR_Obs_SD" = sd(wr_est),

               "Median Observed Ties Proportion" = summary(p_tie_observed)[3],

               "Power"=Power, "Median Lower CI" =summary(lower_ci_bca)[3], "Median Upper CI" = summary(upper_ci_bca)[3],
               "Median CI Width" = summary(CI_width)[3]
               )
  return(mylist)
  WR<-multi_return()


}

