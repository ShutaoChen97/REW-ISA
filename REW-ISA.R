REWISA <- function(FPKM_IP, FPKM_INPUT, Ratio, weight, optimal_thr_row, optimal_thr_col, 
                   optimization, repeat_num, thr_row_interval, row_step,  
                   thr_col_interval, col_step, 
                   fixed_side, side, row_fix_temp, col_fix_temp){
  
  # optimization: Logical variables. If TRUE, turn on threshold optimization
  # repeat_num indicates the number of times to run REW-ISA repeatedly under each pair of threshold parameter settings.
  # thr_row_interval and thr_col_interval represent the selection range of row and column thresholds.
  # row_step and col_step indicate that the row and column threshold is within the step size of the selection.
  # For the setting of the line threshold, it is recommended to change it in the range of 1 to 3 in steps of 0.1.
  # For the setting of the column threshold, it is recommended to change it in steps of 0.05 within 0.1 to 1.5.
  # fixed_side: Usually set to FALSE. Set to TRUE if you want to manually fix the threshold on one side.
  # side: The value 1 indicates a fixed row threshold; the value 2 indicates a fixed column threshold.
  # row_fix_temp: Fixed row threshold value.
  # col_fix_temp: Fixed column threshold value.
  
  # Output of REW-ISA parameter optimization result:
  # ASwC: In each repeated calculation, the Average Similarity within Clusters three-dimensional array calculated for each pair of threshold combinations.
  # SDwC: In each repeated calculation, the Standard Deviation within Clusters three-dimensional array calculated for each pair of threshold combinations.
  # LFB_num: In repeated experiments, a three-dimensional array of LFB numbers generated under each pair of threshold combinations
  # ASwC_mean and SDwC_mean: The average value of each repeated calculation result in each pair of threshold combinations.
  # LFB_num_mode: Under the combination of each pair of thresholds, the mode of the number of LFB is generated.
  # find_TR: Optimized row threshold.
  # find_TC: Optimized col threshold.
  # LFB_number: The optimal number of LFB after optimization.
  
  
  # Import necessary packages
  library(biclust)
  library(isa2)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  
  # Drawing font and color settings
  windowsFonts(RMN = windowsFont("Times New Roman"))
  c_pal <- c("#000066")
  
  if (missing(optimization)) {
    stop("You must specify whether to enable the parameter optimization process.
         Select TRUE or FALSE for optimization.")
  }
  if ((missing(thr_row_interval) && missing(optimal_thr_row) && missing(row_fix_temp))) {
    stop("Missing row threshold related parameters")
  }
  if ((missing(thr_col_interval) && missing(optimal_thr_col) && missing(col_fix_temp))) {
    stop("Missing col threshold related parameters")
  }
  if (optimization == TRUE && fixed_side == TRUE) {
    stop("Conflict between optimized parameters and fixed one-sided")
  }
  
  if ((missing(Ratio) || missing(weight))) {
    FPKM_IP <- FPKM_IP + 0.001
    FPKM_INPUT <- FPKM_INPUT + 0.001
    data_sum <- FPKM_IP + FPKM_INPUT
    Ratio <- FPKM_IP / data_sum
    weight <- log2(data_sum + 1)
  }
  len_row <- nrow(Ratio)
  len_col <- ncol(Ratio)
  
  row_min <- as.numeric(apply(Ratio, 1, min))
  row_max <- as.numeric(apply(Ratio, 1, max))
  row_dis <- row_max - row_min
  row_min_matrix <- matrix(rep(row_min,len_col), ncol=len_col)
  row_dis_matrix <- matrix(rep(row_dis,len_col), ncol=len_col)
  test_1 <- (Ratio - row_min_matrix) / row_dis_matrix
  test_1 <- t(test_1)
  col_min <- apply(Ratio, 2, min)
  col_max <- apply(Ratio, 2, max)
  col_dis <- col_max - col_min
  col_min_matrix <- t(matrix(rep(col_min,len_row), ncol=len_row))
  col_dis_matrix <- t(matrix(rep(col_dis,len_row), ncol=len_row))
  test_2 <- (Ratio - col_min_matrix) / col_dis_matrix
  nm_test <- list()
  nm_test[[1]] <- test_1
  nm_test[[2]] <- test_2
  names(nm_test)[1] <- paste("Er")
  names(nm_test)[2] <- paste("Ec")
  nm_w <- list()
  nm_w[[1]] <- nm_test[[1]] * t(weight)
  nm_w[[2]] <- nm_test[[2]] * weight
  names(nm_w)[1] <- paste("Er")
  names(nm_w)[2] <- paste("Ec")
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  
  if (fixed_side == FALSE) {
    # Grid search threshold optimization
    if (optimization == TRUE) {
      if(repeat_num < 1){
        stop("The number of repetitions must be greater than 0")
      }
      if ((missing(thr_row_interval) || missing(thr_col_interval))) {
        stop("Incorrect interval, please check interval")
      }
      if ((missing(row_step) || missing(col_step))) {
        stop("Incorrect row_step or col_step, please check step!")
      }
      
      thr_row_initial <- thr_row_interval[1]
      thr_row_length <- length(thr_row_interval)
      thr_row_num <- thr_row_length
      thr_row_step <- row_step
      thr_row_final <- thr_row_interval[thr_row_length]
      
      thr_col_initial <- thr_col_interval[1]
      thr_col_length <- length(thr_col_interval)
      thr_col_num <- thr_col_length
      thr_col_step <- col_step
      thr_col_final <- thr_col_interval[thr_col_length]
      
      ASwC <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      SDwC <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      LFB_num <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      ASwC_mean <- array(dim = c(thr_row_num, thr_col_num))
      SDwC_mean <- array(dim = c(thr_row_num, thr_col_num))
      LFB_num_mode <- array(dim = c(thr_row_num, thr_col_num))
      
      for (row_index in 1:thr_row_num) {
        thr_row_temp <- thr_row_interval[row_index]
        cat("\nthr_row_position:", row_index)
        cat("\nthr_row:", thr_row_temp)
        
        for (cycle in 1:repeat_num) {
          cat("\nrepeat:", cycle)
          epoch <- 1
          for (thr in seq(thr_col_initial, thr_col_final, thr_col_step)) {
            cat("\nthr_col:", thr)
            ## Random seeds
            seeds <- generate.seeds(length=len_row, count=100)
            isares <- isa.iterate(nm_w, row.seeds = seeds,
                                  thr.row = thr_row_temp, thr.col = thr)
            ## Eliminate duplicated modules
            isares.unique <- isa.unique(nm_w, isares)
            ## Filter out not robust ones
            isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
            bc <- isa.biclust(isares2)
            # Cluster analysis
            isa_row <- bc@RowxNumber
            isa_col <- bc@NumberxCol
            temp_pcc_vec_col <- vector()
            temp_pcc_post_col <- 1
            mean <- vector()
            std <- vector()
            cluster_number <- ncol(isa_row)
            if(cluster_number < 1) {
              ASwC[cycle, epoch, row_index] = 0
              SDwC[cycle, epoch, row_index] = 0
            }
            if(cluster_number >= 1) {
              for (temp_num in 1:cluster_number) {
                r <- isa_row[, temp_num]
                c <- isa_col[temp_num, ]
                rr <- list()
                i <- 1
                for(var in 1:len_row) {
                  if(r[var]==TRUE) {
                    rr[i] <- var
                    i <- i+1
                  }
                }
                cc <- list()
                i <- 1
                for(var in 1:len_col) {
                  if(c[var]==TRUE) {
                    cc[i] <- var
                    i <- i+1
                  }
                }
                length_c <- length(cc)
                length_r <- length(rr)
                # Calculate ASwC
                for (var_1 in cc) {
                  for (var_2 in cc) {
                    temp_vector_1 <- vector()
                    temp_vector_2 <- vector()
                    temp_post <- 1
                    for (var_3 in rr) {
                      temp_vector_1[temp_post] <- nm_w[[2]][var_3, var_1]
                      temp_vector_2[temp_post] <- nm_w[[2]][var_3, var_2]
                      temp_post <- temp_post + 1
                    }
                    temp_pcc_vec_col[temp_pcc_post_col] <- cor(temp_vector_1, temp_vector_2, method='pearson')
                    temp_pcc_post_col <- temp_pcc_post_col + 1
                  }
                }
                # Calculate SDwC
                temp <- 0
                for (var_1 in rr) {
                  for (var_2 in cc) {
                    temp <- temp + nm_w[[2]][var_1, var_2]
                  }
                }
                mean_temp <- temp / (length_c * length_r)
                mean[temp_num] <- mean_temp
                count <- 0
                for (var_1 in rr) {
                  for (var_2 in cc) {
                    temp <- (nm_w[[2]][var_1, var_2] - mean_temp)^2
                    count <- count + temp
                  }
                }
                std_temp <- sqrt(count / (length_c * length_r))
                std[temp_num] <- std_temp
              }
              ASwC[cycle, epoch, row_index] <- mean(abs(temp_pcc_vec_col))
              SDwC[cycle, epoch, row_index] <- mean(std)
            }
            # Statistics cluster number
            LFB_num[cycle, epoch, row_index] = cluster_number
            epoch <- epoch + 1
          }
        }
      }
      ASwC[is.na(ASwC)] <- 1
      # calculate the mean of ASwC, SDwC and the mode of LFB_num
      for (var_1 in 1:thr_row_num) {
        for (var_2 in 1:thr_col_num) {
          ASwC_mean[var_1, var_2] <- mean(ASwC[, var_2, var_1])
          SDwC_mean[var_1, var_2] <- mean(SDwC[, var_2, var_1])
          LFB_num_mode[var_1, var_2] <- getmode(LFB_num[, var_2, var_1])
        }
      }
      SDAS <- SDwC_mean + ASwC_mean
      # Find out the best LFB_number
      LFB_number_filter <- table(LFB_num_mode)
      LFB_number_filter <- as.data.frame(LFB_number_filter)
      LFB_number_filter <- LFB_number_filter[order(LFB_number_filter[, 2], decreasing = T), ]
      mode_temp <- 1
      for (var_1 in 1:nrow(LFB_number_filter)) {
        if(LFB_number_filter[mode_temp, 1] == 0 ||
           LFB_number_filter[mode_temp, 1] == 1 || 
           LFB_number_filter[mode_temp, 1] == 2) {
          mode_temp <- mode_temp + 1
        }
      }
      LFB_number <- as.numeric(as.character(LFB_number_filter[mode_temp, 1]))
      LFB_logical <- LFB_num_mode == LFB_number
      LFB_select <- LFB_logical * LFB_num_mode
      
      # Establish Sliding window
      LFB_sliding <- matrix(1, nrow = 2, ncol = 2)
      sliding_score <- array(dim = c((thr_row_num - 1), (thr_col_num - 1)))
      sliding_SDAS <- array(dim = c((thr_row_num - 1), (thr_col_num - 1)))
      for (var_1 in 1:(thr_row_num - 1)) {
        for (var_2 in 1:(thr_col_num - 1)) {
          sliding_score[var_1, var_2] <- 
            min(LFB_select[var_1 : (var_1 + 1), var_2 : (var_2 + 1)])
          sliding_SDAS[var_1, var_2] <- 
            mean(SDAS[var_1 : (var_1 + 1), var_2 : (var_2 + 1)])
        }
      }
      sliding_select <- sliding_score * sliding_SDAS
      sliding_select[sliding_select < 0.01] <- 20000
      sliding_row <- which(sliding_select == min(sliding_select), arr.ind = T)[1, 1]
      sliding_col <- which(sliding_select == min(sliding_select), arr.ind = T)[1, 2]
      sliding_row <- as.numeric(sliding_row)
      sliding_col <- as.numeric(sliding_col)
      
      SDAS_min <- SDAS[sliding_row:(sliding_row + 1), sliding_col:(sliding_col + 1)]
      SDAS_row <- which(SDAS_min == min(SDAS_min), arr.ind = T)[1, 1]
      SDAS_col <- which(SDAS_min == min(SDAS_min), arr.ind = T)[1, 2]
      SDAS_row <- as.numeric(SDAS_row)
      SDAS_col <- as.numeric(SDAS_col)
      thr_row_find_index <- sliding_row + SDAS_row - 1
      thr_col_find_index <- sliding_col + SDAS_col - 1
      thr_row_find <- thr_row_interval[thr_row_find_index]
      thr_col_find <- thr_col_interval[thr_col_find_index]
      
      find_TR <<- thr_row_find
      find_TC <<- thr_col_find
      LFB_number <<- LFB_number
      ASwC <<- ASwC
      SDwC <<- SDwC
      LFB_num <<- LFB_num
      ASwC_mean <<- ASwC_mean
      SDwC_mean <<- SDwC_mean
      LFB_num_mode <<- LFB_num_mode
      
      cat("\n\nThe results obtained by REW-ISA are as follows:")
      cat("\noptimal_thr_row:", find_TR)
      cat("\noptimal_thr_col:", find_TC)
      cat("\noptimal_LFB_number:", LFB_number)
    }
    
    
    
    # Finding LFBs under optimized threshold
    if (optimization == FALSE) {
      if ((missing(optimal_thr_row) && missing(optimal_thr_col))) {
        stop("No optimal thresholds, please optimize the row and column thresholds first")
      }
      if ((!missing(optimal_thr_row) && !missing(optimal_thr_col))) {
        # Random seeds
        seeds <- generate.seeds(length = len_row, count = 100)
        isares <- isa.iterate(nm_w, row.seeds = seeds,
                              thr.row = optimal_thr_row, thr.col = optimal_thr_col)
        # Eliminate duplicated modules
        isares.unique <- isa.unique(nm_w, isares)
        # Filter out not robust ones
        isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
        bc <- isa.biclust(isares2)
        LFB <<- bc
      }
    }
  }
  
  
  
  if (fixed_side == TRUE) {
    # draw boxplot
    draw_boxplot <- function(test_matrix, repeat_num,
                             thr_initial, thr_num, thr_step,
                             xlab_name, ylab_name){
      thr_final <- thr_initial + (thr_num - 1) * thr_step
      test_matrix_draw <- as.matrix(test_matrix)
      box_test <- matrix(nrow = repeat_num, ncol = (thr_num + 1))
      box_test[1:repeat_num, 1] <- "REW-ISA"
      box_test <- as.data.frame(box_test)
      box_test[2:(thr_num + 1)] = as.numeric(unlist(box_test[2:(thr_num + 1)]))
      box_test[1:repeat_num, 2:(thr_num + 1)] <- as.numeric(test_matrix_draw)
      colnames(box_test) <- seq(1, (thr_num + 1), 1)
      name <- vector()
      for (var_draw in 1:thr_num) {
        name[var_draw] <- thr_initial + (var_draw - 1) * thr_step
      }
      colnames(box_test) <- c("method", name)
      w_r_score_mean <- apply(test_matrix_draw, 2, mean)
      box_t <- melt(box_test, id.vars = c("method"))
      p1 <- ggplot(box_t) + 
        geom_boxplot(aes(x = variable,y = value,fill = method),width = 0.6,
                     position = position_dodge(0.8),outlier.size = 0,outlier.color = "white")
      p2 <- p1  +
        theme(legend.key = element_blank()) + xlab(xlab_name) + ylab(ylab_name) +
        theme(axis.text.x = element_text(angle = 30,hjust = .5, vjust = .5)) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(plot.title=element_text(hjust=0.5))  + 
        theme(axis.title=element_text(size=12),)
      p3 <- p2 + scale_fill_manual(values = c_pal) +
        scale_color_manual(values = c_pal) + 
        theme(axis.text = element_text(family = "RMN", size = 14), 
              axis.title = element_text(family = "RMN", size = 19)) + 
        theme(legend.title = element_blank(), 
              legend.text = element_text(family = "RMN", size = 21)) +
        scale_x_discrete(breaks = seq(thr_initial, thr_final, 3*thr_step)) + 
        guides(fill = FALSE)
      return(p3)
    }
    
    # draw distribution
    draw_distribution <- function(test_matrix){
      test_matrix_draw <- as.matrix(test_matrix)
      cluster_num_min <- min(test_matrix_draw)
      cluster_num_max <- max(test_matrix_draw)
      a <- as.data.frame(table(test_matrix_draw))
      a2_test <- as.matrix(a)
      a2_test = apply(a2_test, 2, as.numeric)
      aa_row <- nrow(a2_test)
      aa <- matrix(nrow = aa_row, ncol = 3)
      aa[1:nrow(a2_test), 1] <- "REW-ISA"
      aa <- as.data.frame(aa)
      aa[, 2:3] <- as.numeric(unlist(aa[, 2:3]))
      aa[1:nrow(a2_test), 2:3] <- a2_test
      colnames(aa) <- c("method", "variable", "value")
      draw_dis <- ggplot(aa, aes(x = variable, y = value, fill = method)) +
        geom_bar(position = "dodge", stat = "identity") +
        theme(legend.key = element_blank()) + xlab("Number of clusters") +
        ylab("Frequency") + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position = "bottom") + scale_fill_manual(values = c_pal) + 
        scale_color_manual(values = c_pal) + 
        theme(axis.text = element_text(family = "RMN", size = 14), 
              axis.title = element_text(family = "RMN", size = 19)) + 
        theme(legend.title = element_blank(), 
              legend.text = element_text(family = "RMN", size = 21)) +
        guides(fill = FALSE)
      return(draw_dis)
    }
    
    if (side != 1 && side != 2) {
      stop("side must be a numeric value of 1 or 2")
    }
    if (side == 1) {
      if (missing(thr_col_interval) || missing(col_step) || missing(row_fix_temp)) {
        stop("Missing column threshold range or step")
      }
      thr_initial <- thr_col_interval[1]
      thr_length <- length(thr_col_interval)
      thr_num <- thr_length
      thr_step <- col_step
      thr_final <- thr_col_interval[thr_length]
      post <- 1
      cluster_post <- 1
      ASwC_col <- array(dim = c(repeat_num, thr_num))
      SDwC_col <- array(dim = c(repeat_num, thr_num))
      w_col_temp_num_cluster <- array(dim = c(repeat_num, thr_num))
      w_col_cluster_distribution <- vector()
      for (cycle in 1:repeat_num) {
        epoch <- 1
        for (thr in seq(thr_initial, thr_final, thr_step)) {
          ## Random seeds
          seeds <- generate.seeds(length = len_row_final, count = 100)
          isares <- isa.iterate(nm_w, row.seeds = seeds, thr.row = row_fix_temp, thr.col = thr)
          ## Eliminate duplicated modules
          isares.unique <- isa.unique(nm_w, isares)
          ## Filter out not robust ones
          isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
          bc <- isa.biclust(isares2)
          # Cluster analysis
          isa_row <- bc@RowxNumber
          isa_col <- bc@NumberxCol
          temp_pcc_vec_col <- vector()
          temp_pcc_post_col <- 1
          mean <- vector()
          std <- vector()
          cluster_number <- ncol(isa_row)
          w_col_cluster_distribution[cluster_post] <- cluster_number
          cluster_post <- cluster_post + 1
          if(cluster_number < 1) {
            ASwC_col[post, epoch] = 0
            SDwC_col[post, epoch] = 0
          }
          if(cluster_number >= 1) {
            for (temp_num in 1:cluster_number) {
              r <- isa_row[, temp_num]
              c <- isa_col[temp_num, ]
              rr <- list()
              i <- 1
              for(var in 1:len_row_final) {
                if(r[var]==TRUE) {
                  rr[i] <- var
                  i <- i+1
                }
              }
              cc <- list()
              i <- 1
              for(var in 1:len_col_final) {
                if(c[var]==TRUE) {
                  cc[i] <- var
                  i <- i+1
                }
              }
              length_c <- length(cc)
              length_r <- length(rr)
              # Calculate ASwC
              for (var_1 in cc) {
                for (var_2 in cc) {
                  temp_vector_1 <- vector()
                  temp_vector_2 <- vector()
                  temp_post <- 1
                  for (var_3 in rr) {
                    temp_vector_1[temp_post] <- nm_w[[2]][var_3, var_1]
                    temp_vector_2[temp_post] <- nm_w[[2]][var_3, var_2]
                    temp_post <- temp_post + 1
                  }
                  temp_pcc_vec_col[temp_pcc_post_col] <- cor(temp_vector_1, temp_vector_2, method='pearson')
                  temp_pcc_post_col <- temp_pcc_post_col + 1
                }
              }
              # Calculate SDwC
              temp <- 0
              for (var_1 in rr) {
                for (var_2 in cc) {
                  temp <- temp + nm_w[[2]][var_1, var_2]
                }
              }
              mean_temp <- temp / (length_c * length_r)
              mean[temp_num] <- mean_temp
              count <- 0
              for (var_1 in rr) {
                for (var_2 in cc) {
                  temp <- (nm_w[[2]][var_1, var_2] - mean_temp)^2
                  count <- count + temp
                }
              }
              std_temp <- sqrt(count / (length_c * length_r))
              std[temp_num] <- std_temp
            }
            ASwC_col[post, epoch] <- mean(abs(temp_pcc_vec_col))
            SDwC_col[post, epoch] = mean(std)
          }
          # Statistics cluster number
          w_col_temp_num_cluster[post, epoch] = cluster_number
          epoch <- epoch + 1
        }
        post <- post + 1
      }
      
      # Visualization
      # draw SDwC_col
      p_col_std <- draw_boxplot(SDwC_col, repeat_num = repeat_num,
                                thr_initial = thr_initial,
                                thr_num = thr_num, thr_step = thr_step,
                                xlab_name = expression(italic(T[C])),
                                ylab_name = "SDwC")
      
      # draw ASwC_col
      p_col_pearson <- draw_boxplot(ASwC_col, repeat_num = repeat_num,
                                    thr_initial = thr_initial,
                                    thr_num = thr_num, thr_step = thr_step,
                                    xlab_name = expression(italic(T[C])),
                                    ylab_name = "ASwC")
      
      # draw w_col_temp_num_cluster
      p_col_number <- draw_boxplot(w_col_temp_num_cluster,
                                   repeat_num = repeat_num,
                                   thr_initial = thr_initial,
                                   thr_num = thr_num, thr_step = thr_step,
                                   xlab_name = expression(italic(T[C])),
                                   ylab_name = "Cluster number")
      
      # draw w_col_cluster_distribution
      cluster_dis_col <- draw_distribution(w_col_cluster_distribution)
      
      # draw total
      draw_col <- ggarrange(cluster_dis_col, p_col_number,
                            p_col_std, p_col_pearson,
                            labels = c("A", "B", "C", "D"),
                            common.legend = FALSE, legend = "bottom")
      print(draw_col)
    }
    
    
    if (side == 2) {
      if (missing(thr_row_interval) || missing(row_step) || missing(col_fix_temp)) {
        stop("Missing row threshold range or step")
      }
      thr_initial <- thr_row_interval[1]
      thr_length <- length(thr_row_interval)
      thr_num <- thr_length
      thr_step <- row_step
      thr_final <- thr_row_interval[thr_length]
      
      post <- 1
      cluster_post <- 1
      ASwC_row <- array(dim=c(repeat_num, thr_num))
      SDwC_row <- array(dim=c(repeat_num, thr_num))
      w_row_temp_num_cluster <- array(dim=c(repeat_num, thr_num))
      w_row_cluster_distribution <- vector()
      for (cycle in 1:repeat_num) {
        epoch <- 1
        for (thr in seq(thr_initial, thr_final, thr_step)) {
          ## Random seeds
          seeds <- generate.seeds(length=len_row_final, count=100)
          isares <- isa.iterate(nm_w, row.seeds = seeds, thr.row = thr, thr.col = col_fix_temp)
          ## Eliminate duplicated modules
          isares.unique <- isa.unique(nm_w, isares)
          ## Filter out not robust ones
          isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
          bc <- isa.biclust(isares2)
          # Cluster analysis
          isa_row <- bc@RowxNumber
          isa_col <- bc@NumberxCol
          temp_pcc_vec_col <- vector()
          temp_pcc_post_col <- 1
          mean <- vector()
          std <- vector()
          cluster_number <- ncol(isa_row)
          w_row_cluster_distribution[cluster_post] <- cluster_number
          cluster_post <- cluster_post + 1
          if(cluster_number < 1) {
            ASwC_row[post, epoch] = 0
            SDwC_row[post, epoch] = 0
          }
          if(cluster_number >= 1) {
            for (temp_num in 1:cluster_number) {
              r <- isa_row[, temp_num]
              c <- isa_col[temp_num, ]
              rr <- list()
              i <- 1
              for(var in 1:len_row_final) {
                if(r[var]==TRUE) {
                  rr[i] <- var
                  i <- i+1
                }
              }
              cc <- list()
              i <- 1
              for(var in 1:len_col_final) {
                if(c[var]==TRUE) {
                  cc[i] <- var
                  i <- i+1
                }
              }
              length_c <- length(cc)
              length_r <- length(rr)
              # Calculate ASwC
              for (var_1 in cc) {
                for (var_2 in cc) {
                  temp_vector_1 <- vector()
                  temp_vector_2 <- vector()
                  temp_post <- 1
                  for (var_3 in rr) {
                    temp_vector_1[temp_post] <- nm_w[[2]][var_3, var_1]
                    temp_vector_2[temp_post] <- nm_w[[2]][var_3, var_2]
                    temp_post <- temp_post + 1
                  }
                  temp_pcc_vec_col[temp_pcc_post_col] <- cor(temp_vector_1, temp_vector_2, method='pearson')
                  temp_pcc_post_col <- temp_pcc_post_col + 1
                }
              }
              # Calculate SDwC
              temp <- 0
              for (var_1 in rr) {
                for (var_2 in cc) {
                  temp <- temp + nm_w[[2]][var_1, var_2]
                }
              }
              mean_temp <- temp / (length_c * length_r)
              mean[temp_num] <- mean_temp
              count <- 0
              for (var_1 in rr) {
                for (var_2 in cc) {
                  temp <- (nm_w[[2]][var_1, var_2] - mean_temp)^2
                  count <- count + temp
                }
              }
              std_temp <- sqrt(count / (length_c * length_r))
              std[temp_num] <- std_temp
            }
            ASwC_row[post, epoch] <- mean(abs(temp_pcc_vec_col))
            SDwC_row[post, epoch] = mean(std)
          }
          # Statistics cluster number
          w_row_temp_num_cluster[post, epoch] = cluster_number
          epoch <- epoch + 1
        }
        post <- post + 1
      }
      
      # Visualization
      # draw SDwC_row
      p_row_std <- draw_boxplot(SDwC_row, repeat_num = repeat_num,
                                thr_initial = thr_initial,
                                thr_num = thr_num, thr_step = thr_step,
                                xlab_name = expression(italic(T[R])),
                                ylab_name = "SDwC")
      
      # draw ASwC_row
      p_row_pearson <- draw_boxplot(ASwC_row, repeat_num = repeat_num,
                                    thr_initial = thr_initial,
                                    thr_num = thr_num, thr_step = thr_step,
                                    xlab_name = expression(italic(T[R])),
                                    ylab_name = "ASwC")
      
      # draw w_row_temp_num_cluster
      p_row_number <- draw_boxplot(w_row_temp_num_cluster,
                                   repeat_num = repeat_num,
                                   thr_initial = thr_initial,
                                   thr_num = thr_num, thr_step = thr_step,
                                   xlab_name = expression(italic(T[R])),
                                   ylab_name = "Cluster number")
      
      # draw w_row_cluster_distribution
      cluster_dis_row <- draw_distribution(w_row_cluster_distribution)
      
      # draw total
      draw_row <- ggarrange(cluster_dis_row, p_row_number,
                            p_row_std, p_row_pearson,
                            labels = c("A", "B", "C", "D"),
                            common.legend = FALSE, legend = "bottom")
      print(draw_row)
    }
  }
}

