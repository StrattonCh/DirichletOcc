label.switch.sfw <- function(x, dp.con = "alpha", dp.latent = "z", cluster.on = "theta",
                              k.fixed = "auto"){
  # helpers
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  getcols <- function(y, param) y[,grepl(param, colnames(y))]
  
  # identify parameters
  params <- unique(sapply(strsplit(colnames(x), split = "[[]"), function(x) x[[1]][1]))
  params_nondp <- params[-which(params %in% c(dp.con, dp.latent))]
  
  # isolate each parameter
  params.list <- list()
  for(i in 1:length(params)){
    params.list[[params[i]]] <- getcols(x, params[i])
  }
  
  # step 1
  k0 <- apply(params.list[[dp.latent]], 1, function(x) length(unique(x)))
  
  # step 2
  khat <- getmode(k0)
  if(k.fixed != "auto"){
    khat <- k.fixed
  }
  
  # step 3
  ## remove iterations where k0 != khat
  x_k0 <- x[which(k0 == khat),]
  params.list <- list()
  for(i in 1:length(params)){
    params.list[[params[i]]] <- getcols(x_k0, params[i])
  }
  
  ## relabel
  params.list.reorder <- list()
  for(i in 1:length(params_nondp)) params.list.reorder[[i]] <- matrix(0, nrow = nrow(x_k0), ncol = khat)
  names(params.list.reorder) <- params_nondp
  latent.reorder <- params.list[[dp.latent]]
  for(i in 1:nrow(x_k0)){
    latent_i <- params.list[[dp.latent]][i,]
    latent_i_uniq <- unique(latent_i)
    for(j in 1:length(params_nondp)){
      params.list.reorder[[params_nondp[j]]][i,] <- params.list[[params_nondp[j]]][i,latent_i_uniq]
    }
    latent.reorder[i,] <- as.numeric(factor(latent_i, levels = latent_i_uniq))
  }
  
  # step 4
  data_matrix <- sapply(cluster.on, function(x){
    c(t(as.matrix(params.list.reorder[[x]], nrow = nrow(x_k0))))
  })
  kmeans <- kmeans(data_matrix, centers = khat, nstart = 25)
  rho_mat <- matrix(kmeans$cluster, ncol = khat, byrow = TRUE)
  
  # step 5
  isperm <- apply(rho_mat, 1, function(x) all(1:khat %in% x))
  if(all(!isperm)) stop("No classification vectors are a permuation of 1:K. Run MCMC for more iterations.")
  rho_mat <- rho_mat[isperm,]
  latent.reorder <- latent.reorder[isperm,]
  for(i in 1:length(params.list.reorder)) params.list.reorder[[i]] <- params.list.reorder[[i]][isperm,]
  con <- as.matrix(x_k0[isperm,dp.con], ncol = 1);colnames(con) <- dp.con
  
  # step 6
  for(i in 1:sum(isperm)){
    for(j in 1:length(params.list.reorder)){
      params.list.reorder[[j]][i,] <- params.list.reorder[[j]][i,order(rho_mat[i,])]
    }
    
    tmp <- latent.reorder[i,]
    for(k in 1:khat){
      tmp[which(latent.reorder[i,] == k)] <- rho_mat[i, k]
    }
    latent.reorder[i,] <- tmp
  }
  
  head(rho_mat)
  head(params.list.reorder[[2]])
  head(latent.reorder[,1:6])
  head(data$y)
  
  for(j in 1:length(params.list.reorder)){
    colnames(params.list.reorder[[j]]) <- paste0(names(params.list.reorder)[j], "[", 1:khat, "]")
  }
  
  # out 
  out <- cbind(
    con,
    do.call("cbind", params.list.reorder), 
    latent.reorder
  )
  
  return(out)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dotplot <- function(x, groups, prob, bins = 30, mean = NULL, sigma2 = NULL, size = 3){
  hist_obj <- hist(x, plot = FALSE, breaks = 30)
  break_mat <- cbind(
    hist_obj$breaks[1:(length(hist_obj$breaks)-1)],
    hist_obj$breaks[2:(length(hist_obj$breaks))]
  )
  
  isbetween <- function(x, mat) {
    x <- matrix(rep(x, nrow(mat)), ncol = 1)
    tmp <- cbind(x, mat)
    which(apply(tmp, 1, function(x) x[1] > x[2] & x[1] < x[3]))
  }
  
  bin <- sapply(x, function(y) isbetween(y, break_mat))
  
  tmp <- tibble(
    data = x, 
    groups = groups, 
    prob = prob,
    bin = bin
  ) %>% 
    mutate(x = hist_obj$mids[bin]) %>%
    # dplyr::arrange(bin) %>%
    group_by(bin) %>%
    mutate(y = 1:n()) %>%
    ungroup() %>%
    mutate(prob = I(prob))
  
  p <- ggplot(tmp) +
    geom_point(aes(x = x, y = y, col = groups, alpha = prob), size = size) + 
    theme_bw()
  
  if(!is.null(mean)){
    binwidth <- unique(diff(hist_obj$breaks))
    colors <- gg_color_hue(length(unique(groups)))
    for(i in 1:length(mean)){
      # p <- p + 
      # stat_function(fun = function(x) dnorm(x, mean[i], sqrt(sigma2[i])) * table(groups)[i] * binwidth, color = colors[i])
      
      df <- data.frame(
        x = rnorm(1000, mean[i], sqrt(sigma2[i]))
      ) %>%
        mutate(
          y = dnorm(x, mean[i], sqrt(sigma2[i])) * binwidth * table(groups)[i]
        )
      p <- p + 
        geom_line(data = df, aes(x = x, y = y), color = colors[i])
    }
  }
  
  return(p)
  
}
