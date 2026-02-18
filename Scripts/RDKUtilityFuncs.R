# Utility -----------------------------------------------------------------

data_clean_rdk <- function(x){

  if('Participant Public ID' %in% colnames(x)){
    y <- x %>%
      as.data.frame() %>%
      rename(ID = `Participant Public ID`,
             Trial = `Trial number`,
             RT = `Total time for reference selection`,
             mistake = 'confidenceMistake',
             conf = 'confidenceSelection',
             type = `Trial Type`) %>%
      dplyr::select(Trial, ID, correct, RT, coherence, conf, dotDirection, referenceSelection, mistake) %>%
      filter(Trial != is.na(Trial),
             Trial > 19) %>%
      group_by(ID) %>%
      mutate(Trial = 1:320,
             type = c(rep('calibration', 120), rep('main', 200))) %>%
      group_by(ID, type) %>%
      mutate( kmed = ifelse(coherence == max(coherence), 'kmed x 2',
                            ifelse(coherence == min(coherence), 'kmed x 0.5', 'kmed')),
              kmed = factor(kmed, levels = c('kmed x 0.5', 'kmed', 'kmed x 2'), ordered = T),
              corAdjusted = ifelse(dotDirection == 3.142 & referenceSelection == 1, 1,
                                   ifelse(dotDirection == 0.000 & referenceSelection == 2, 1, 0))
      )
  } else {

    y <- x %>%
      group_by(ID) %>%
      as.data.frame() %>%
      rename(
        Trial = `trialNumber`,
        RT = referenceTotalTime,
        mistake = 'confidenceMistake',
        conf = 'confidenceSelection',
        type = name) %>%
      dplyr::select(Trial, ID, RT, type, correct, coherence, conf, dotDirection, referenceSelection, mistake) %>%
      filter(Trial != is.na(Trial)) %>%
      unique() %>%
      group_by(ID) %>%
      filter(type %in% c('calibration', 'main')) %>%
      group_by(ID, type) %>%
      mutate(Trial = 1:length(type),
             kmed = ifelse(coherence == max(coherence), 'kmed x 2',
                           ifelse(coherence == min(coherence), 'kmed x 0.5', 'kmed')),
             kmed = factor(kmed, levels = c('kmed x 0.5', 'kmed', 'kmed x 2'), ordered = T),
             corAdjusted = ifelse(dotDirection == 3.142 & referenceSelection == 1, 1,
                                  ifelse(dotDirection == 0.000 & referenceSelection == 2, 1, 0))
      ) %>%
      as.data.frame()
  }

  y$type <- replace(y$type, grepl('calibration', y$type, fixed=TRUE), 'calibration')

  return(y)
}

check_coherence <- function(x, groupname, legend=T) {

  x <- x %>% filter(group == groupname)

  if(nrow(subset(x %>%
                 dplyr::filter(type == "calibration") %>%
                 dplyr::select(ID, Trial, coherence),
                 coherence > 0.5)) > 1){
    Switch <- T
  } else {
    Switch <- F
  }

  y <- ggplot(x %>%
                dplyr::filter(type == "calibration") %>%
                dplyr::select(ID, Trial, coherence)) +
    geom_line(aes(x = Trial, y = coherence, colour = ID, group = ID)) +
    geom_hline(yintercept = 0.50, colour = "black", linetype="dashed") +
    geom_hline(yintercept = 0.12, colour = "black", linetype="dashed") +
    {if(Switch!=F)geom_text(data = subset(x %>%
                                            dplyr::filter(type == "calibration") %>%
                                            dplyr::select(ID, Trial, coherence),
                                          coherence > 0.5),
                            aes(y = 0.6, x = 100, label = ID),
                            check_overlap = T, fontface = 'bold')}+
    {if(legend!=T)guides(colour = 'none')}+
    labs(x = "Trial number", y = "Coherence", title = groupname) +
    coord_cartesian(ylim=c(0.0, 0.6), xlim = c(0, 120)) +
    theme_bw()+
    theme(text = element_text(size=20))
  return(y)
}

HMetaPrep <- function(x, DD, ref, desc){

  y <- x %>%
    filter(type == 'main',
           conf != 0) %>%
    mutate(conf = factor(conf, levels = c(50, 60, 70, 80, 90, 100), order = T),
           dotDirection = ifelse(dotDirection == 0, '2', '1')) %>%
    group_by(ID, conf, .drop = FALSE) %>%
    filter(dotDirection == DD,
           referenceSelection == ref) %>%
    summarise(nR_S = n())

  if (desc == T){
    y <- y %>% arrange(ID, wt=desc(conf))
  }else{
    y <- y %>% arrange(ID, wt=conf)
  }
  return(y)
}

HMetaGroupPrep <- function(x, groupName) {

  nR_S1raw1 <- HMetaPrep(x %>% filter(group == groupName), '1', 1, T)
  nR_S1raw2 <- HMetaPrep(x %>% filter(group == groupName), '1', 2, F)
  nR_S1raw  <- rbind(nR_S1raw1, nR_S1raw2) %>% arrange(ID)
  nR_S1raw$groupID <-  nR_S1raw %>% group_indices(ID, .groups = T)

  nR_S1  <- matrix(NA, nrow = 12, ncol = length(unique(nR_S1raw$groupID)))
  for (i in 1:length(unique(nR_S1raw$groupID))){
    y <- nR_S1raw  %>%
      filter(groupID == unique(nR_S1raw$groupID)[i]) %>%
      ungroup() %>%
      dplyr::select(nR_S) %>%
      as.data.frame() %>%
      as.vector() %>%
      unlist()
    nR_S1 [,i] <- y
  }

  nR_S2raw1 <- HMetaPrep(x %>% filter(group == groupName), '2', 1, T)
  nR_S2raw2 <- HMetaPrep(x %>% filter(group == groupName), '2', 2, F)
  nR_S2raw <- rbind(nR_S2raw1, nR_S2raw2) %>% arrange(ID)
  nR_S2raw$groupID <-  nR_S2raw  %>% group_indices(ID, .groups = T)

  nR_S2  <- matrix(NA, nrow = 12, ncol = length(unique(nR_S2raw$groupID)))
  for (i in 1:length(unique(nR_S2raw$groupID))){
    y <- nR_S2raw  %>%
      filter(groupID == unique(nR_S2raw$groupID)[i]) %>%
      ungroup() %>%
      dplyr::select(nR_S) %>%
      as.data.frame() %>%
      as.vector()%>%
      unlist()
    nR_S2 [,i] <- y
  }

  return(list(S1mat = nR_S1 %>% as.data.frame(),
              S2mat = nR_S2 %>% as.data.frame()))

}

HMetaPrep_vr <- function(x, stim, ref, desc){

  y <- x %>%
        filter(TrialType == 'Main') %>%
        group_by(PID, ReportedVeryConfident, .drop = FALSE) %>%
        filter(ifelse(DotDirection == 'Up', 1, 0) == stim, ReportedUpwardMovement == ref) %>%
        summarise(nR_S = n()) %>%
        {summarized_data <- .
         bind_rows(
           summarized_data,
           expand_grid(PID = unique(summarized_data$PID), ReportedVeryConfident = c(0, 1)) %>%
           anti_join(summarized_data, by = c("PID", "ReportedVeryConfident")) %>%
           mutate(nR_S = 0)
         )} %>%
         arrange(PID, ReportedVeryConfident)

  if (desc == T)
  {
    y <- y %>% arrange(PID, wt=desc(ReportedVeryConfident))
  }
  else
  {
    y <- y %>% arrange(PID, wt=ReportedVeryConfident)
  }

  return(y)
}

HMetaGroupPrep_vr <- function(x, groupName) {

    nR_S1raw1 <- HMetaPrep_vr(x %>% filter(Group == groupName), 1, 1, T)
    nR_S1raw2 <- HMetaPrep_vr(x %>% filter(Group == groupName), 1, 0, F)
    nR_S1raw  <- rbind(nR_S1raw1, nR_S1raw2) %>% arrange(PID)
    nR_S1raw$groupID <-  nR_S1raw %>% group_indices(PID, .groups = T)
    nR_S1raw  <- nR_S1raw %>% group_by(PID) %>% filter(!n()<4)

    nR_S1  <- matrix(NA, nrow = 4, ncol = length(unique(nR_S1raw$groupID)))
    for (i in 1:length(unique(nR_S1raw$groupID)))
    {
      y <- nR_S1raw  %>%
        filter(groupID == unique(nR_S1raw$groupID)[i]) %>%
        ungroup() %>%
        dplyr::select(nR_S) %>%
        as.data.frame() %>%
        as.vector() %>%
        unlist()
      nR_S1 [,i] <- y
    }

    nR_S2raw1 <- HMetaPrep_vr(x %>% filter(Group == groupName), 0, 1, T)
    nR_S2raw2 <- HMetaPrep_vr(x %>% filter(Group == groupName), 0, 0, F)
    nR_S2raw <- rbind(nR_S2raw1, nR_S2raw2) %>% arrange(PID)
    nR_S2raw$groupID <-  nR_S2raw  %>% group_indices(PID, .groups = T)
    nR_S2raw  <- nR_S2raw %>% group_by(PID) %>% filter(!n()<4)

    nR_S2  <- matrix(NA, nrow = 4, ncol = length(unique(nR_S2raw$groupID)))
    for (i in 1:length(unique(nR_S2raw$groupID))){
      y <- nR_S2raw  %>%
        filter(groupID == unique(nR_S2raw$groupID)[i]) %>%
        ungroup() %>%
        dplyr::select(nR_S) %>%
        as.data.frame() %>%
        as.vector()%>%
        unlist()
      nR_S2 [,i] <- y
    }

  return(list(S1mat = nR_S1 %>% as.data.frame(),
              S2mat = nR_S2 %>% as.data.frame()))
}


HMeta_post_clean <- function(x, groupname, trace = T){

  Value <- summary(x)
  stat  <- data.frame(mean = Value$statistics[,"Mean"])
  stat %<>%
    rownames_to_column(var = "name")

  Value      <- gelman.diag(x, confidence = 0.95)
  Rhat       <- data.frame(conv = Value$psrf)

  hdi <- data.frame(HPDinterval(x, prob = 0.95))
  hdi %<>%
    rownames_to_column(var = "name")

  fit <- stat %>%
    cbind(lower = hdi$lower,
          upper = hdi$upper,
          Rhat = Rhat[,1])

  mcmc.sample <- ggs(x) %>% mutate(group = groupname)

  if(trace == T){
    traceplot(x)
  }

  return(list(MCMC = mcmc.sample,
              Value = Value,
              Fit  = fit))

}

HMeta_post_clean_meta_d <- function(x, groupname, trace = F){

  Value <- summary(x)
  stat  <- data.frame(mean = Value$statistics[,"Mean"])
  stat %<>%
    rownames_to_column(var = "name")

  hdi <- data.frame(HPDinterval(x, prob = 0.95))
  hdi %<>%
    rownames_to_column(var = "name")

  fit <- stat %>%
    cbind(lower = hdi$lower,
          upper = hdi$upper)

  mcmc.sample <- ggs(x) %>% mutate(group = groupname)

  if(trace == T){
    traceplot(x)
  }

  return(list(MCMC = mcmc.sample,
              Value = Value,
              Fit  = fit))

}


permute_metad_group <- function(df1, df2, cores = 4, nreps = 500){

  doParallel::registerDoParallel(cores = cores)

  permuted_metad <- foreach(i=1:nreps, .combine = rbind) %dopar% {

    r.random   <- matrix(NA, 3, 3)

    nid_of_df1 <- length(unique(df1$ID))
    permuteddf <- rbind(df1, df2)
    IDs        <- unique(permuteddf$ID)
    IDstoUse   <- sample(IDs, nid_of_df1)
    DF1permute <- permuteddf %>% filter(ID %in% IDstoUse) %>% mutate(group = 'DF1')
    DF2permute <- permuteddf %>% filter(!ID %in% IDstoUse) %>% mutate(group = 'DF2')

    DF1run <- HMetaGroupPrep(DF1permute, 'DF1')
    DF2run <- HMetaGroupPrep(DF2permute, 'DF2')

    nR_S1_pDF1<-list(as.data.frame(DF1run[[1]]))
    nR_S2_pDF1<-list(as.data.frame(DF1run[[2]]))
    nR_S1_pDF2<-list(as.data.frame(DF2run[[1]]))
    nR_S2_pDF2<-list(as.data.frame(DF2run[[2]]))

    #hierarchical fit
    outputPermuteDF1 <- try(metad_only_group(nR_S1 = nR_S1_pDF1, nR_S2 = nR_S2_pDF1))
    outputPermuteDF2 <- try(metad_only_group(nR_S1 = nR_S1_pDF2, nR_S2 = nR_S2_pDF2))

    if(!inherits(outputPermuteDF1, 'try-error') &
       !inherits(outputPermuteDF2, 'try-error')) {

      ValueDF1p <- summary(outputPermuteDF1)
      statDF1p <- data.frame(mean = ValueDF1p$statistics[,"Mean"])
      statDF1p %<>%
        rownames_to_column(var = "name")
      ValueDF2p <- summary(outputPermuteDF2)
      statDF2p <- data.frame(mean = ValueDF2p$statistics[,"Mean"])
      statDF2p %<>%
        rownames_to_column(var = "name")

      r.random[1,1] <- statDF1p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[2,1] <- statDF2p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[3,1] <- r.random[2, 1] - r.random[1, 1]

      r.random[1,2] <- statDF1p %>%
        filter(name == 'mu_d1') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[2,2] <- statDF2p %>%
        filter(name == 'mu_d1') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[3,2] <- r.random[2, 2] - r.random[1, 2]

      r.random[1,3] <- statDF1p %>%
        filter(name == 'mu_meta_d') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[2,3] <- statDF2p %>%
        filter(name == 'mu_meta_d') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[3,3] <- r.random[2, 3] - r.random[1, 3]

    } else {

      r.random[1:3, 1:3] <- NA

    }

    data.frame(
      rep = i,
      muGroup1     = r.random[1,1],
      muGroup2     = r.random[2,1],
      mu_diff      = r.random[3,1],
      d1Group1     = r.random[1,2],
      d1Group2     = r.random[2,2],
      d1_diff      = r.random[3,2],
      meta_dGroup1 = r.random[1,3],
      meta_dGroup2 = r.random[2,3],
      meta_d_diff  = r.random[3,3]

    )

  }

  return(permuted_metad)

}

permute_metad_group_mrat <- function(df1, df2, cores = 4, nreps = 500){

  doParallel::registerDoParallel(cores = cores)

  permuted_metad <- foreach(i=1:nreps, .combine = rbind) %dopar% {

    r.random   <- rep(NA, 3)

    nid_of_df1 <- length(unique(df1$ID))
    permuteddf <- rbind(df1, df2)
    IDs        <- unique(permuteddf$ID)
    IDstoUse   <- sample(IDs, nid_of_df1)
    DF1permute <- permuteddf %>% filter(ID %in% IDstoUse) %>% mutate(group = 'DF1')
    DF2permute <- permuteddf %>% filter(!ID %in% IDstoUse) %>% mutate(group = 'DF2')

    DF1run <- HMetaGroupPrep(DF1permute, 'DF1')
    DF2run <- HMetaGroupPrep(DF2permute, 'DF2')

    nR_S1_pDF1<-list(as.data.frame(DF1run[[1]]))
    nR_S2_pDF1<-list(as.data.frame(DF1run[[2]]))
    nR_S1_pDF2<-list(as.data.frame(DF2run[[1]]))
    nR_S2_pDF2<-list(as.data.frame(DF2run[[2]]))

    #hierarchical fit
    outputPermuteDF1 <- try(metad_group(nR_S1 = nR_S1_pDF1, nR_S2 = nR_S2_pDF1))
    outputPermuteDF2 <- try(metad_group(nR_S1 = nR_S1_pDF2, nR_S2 = nR_S2_pDF2))

    if(!inherits(outputPermuteDF1, 'try-error') &
       !inherits(outputPermuteDF2, 'try-error')) {

      ValueDF1p <- summary(outputPermuteDF1)
      statDF1p <- data.frame(mean = ValueDF1p$statistics[,"Mean"])
      statDF1p %<>%
        rownames_to_column(var = "name")
      ValueDF2p <- summary(outputPermuteDF2)
      statDF2p <- data.frame(mean = ValueDF2p$statistics[,"Mean"])
      statDF2p %<>%
        rownames_to_column(var = "name")

      r.random[1] <- statDF1p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[2] <- statDF2p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[3] <- r.random[2] - r.random[1]

    } else {

      r.random[1:3] <- NA

    }

    data.frame(
      rep = i,
      muGroup1     = r.random[1],
      muGroup2     = r.random[2],
      mu_diff      = r.random[3]

    )

  }

  return(permuted_metad)

}

permute_metad_group_vr <- function(df1, df2, cores = 4, nreps = 500){

  doParallel::registerDoParallel(cores = cores)

  permuted_metad <- foreach(i=1:nreps, .combine = rbind) %dopar% {

    r.random   <- matrix(NA, 3, 3)

    nid_of_df1 <- length(unique(df1$PID))
    permuteddf <- rbind(df1, df2)
    IDs        <- unique(permuteddf$PID)
    IDstoUse   <- sample(IDs, nid_of_df1)
    DF1permute <- permuteddf %>% filter(PID %in% IDstoUse) %>% mutate(Group = 'DF1')
    DF2permute <- permuteddf %>% filter(!PID %in% IDstoUse) %>% mutate(Group = 'DF2')

    DF1run <- HMetaGroupPrep_vr(DF1permute, 'DF1')
    DF2run <- HMetaGroupPrep_vr(DF2permute, 'DF2')

    nR_S1_pDF1<-list(as.data.frame(DF1run[[1]]))
    nR_S2_pDF1<-list(as.data.frame(DF1run[[2]]))
    nR_S1_pDF2<-list(as.data.frame(DF2run[[1]]))
    nR_S2_pDF2<-list(as.data.frame(DF2run[[2]]))

    #hierarchical fit
    outputPermuteDF1 <- try(metad_only_group(nR_S1 = nR_S1_pDF1, nR_S2 = nR_S2_pDF1))
    outputPermuteDF2 <- try(metad_only_group(nR_S1 = nR_S1_pDF2, nR_S2 = nR_S2_pDF2))

    if(!inherits(outputPermuteDF1, 'try-error') &
       !inherits(outputPermuteDF2, 'try-error')) {

      ValueDF1p <- summary(outputPermuteDF1)
      statDF1p <- data.frame(mean = ValueDF1p$statistics[,"Mean"])
      statDF1p %<>%
        rownames_to_column(var = "name")
      ValueDF2p <- summary(outputPermuteDF2)
      statDF2p <- data.frame(mean = ValueDF2p$statistics[,"Mean"])
      statDF2p %<>%
        rownames_to_column(var = "name")

      r.random[1,1] <- statDF1p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[2,1] <- statDF2p %>%
        filter(name == 'mu_logMratio') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        exp()
      r.random[3,1] <- r.random[2, 1] - r.random[1, 1]

      r.random[1,2] <- statDF1p %>%
        filter(name == 'mu_d1') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[2,2] <- statDF2p %>%
        filter(name == 'mu_d1') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[3,2] <- r.random[2, 2] - r.random[1, 2]

      r.random[1,3] <- statDF1p %>%
        filter(name == 'mu_meta_d') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[2,3] <- statDF2p %>%
        filter(name == 'mu_meta_d') %>%
        dplyr::select(mean) %>%
        as.numeric() %>%
        abs()
      r.random[3,3] <- r.random[2, 3] - r.random[1, 3]

    } else {

      r.random[1:3, 1:3] <- NA

    }

    data.frame(
      rep = i,
      muGroup1     = r.random[1,1],
      muGroup2     = r.random[2,1],
      mu_diff      = r.random[3,1],
      d1Group1     = r.random[1,2],
      d1Group2     = r.random[2,2],
      d1_diff      = r.random[3,2],
      meta_dGroup1 = r.random[1,3],
      meta_dGroup2 = r.random[2,3],
      meta_d_diff  = r.random[3,3]

    )

  }

  return(permuted_metad)

}

permute_metad_group_vr_within <- function(df1, df2, cores = 4, nreps = 500) {
  # Ensure pairing exists
  common_ids <- intersect(unique(df1$PID), unique(df2$PID))
  if (length(common_ids) == 0L) {
    stop("No overlapping PIDs between df1 and df2; paired permutation is not possible.")
  }

  # Keep only paired subjects
  df1p <- df1[df1$PID %in% common_ids, , drop = FALSE]
  df2p <- df2[df2$PID %in% common_ids, , drop = FALSE]

  # Tag original group and stack
  df1p$Group <- "DF1"
  df2p$Group <- "DF2"
  all_df <- rbind(df1p, df2p)

  # Parallel setup (reproducible)
  doParallel::registerDoParallel(cores = cores)
  # Optional: set.seed() before doRNG to control the whole stream
  if (!requireNamespace("doRNG", quietly = TRUE)) {
    stop("Package 'doRNG' is required for reproducible parallel RNG. Please install it.")
  }
  doRNG::registerDoRNG()  # reproducible %dopar% RNG

  # Main loop
  permuted_metad <- foreach::foreach(
    i = 1:nreps,
    .combine   = rbind,
    .inorder   = FALSE,
    .packages  = c("dplyr")
  ) %dopar% {
    # 50/50 within-subject label swap
    swap <- stats::rbinom(length(common_ids), 1, 0.5) == 1
    swap_map <- data.frame(PID = common_ids, swap = swap)

    permuted <- dplyr::left_join(all_df, swap_map, by = "PID")
    permuted$Group_perm <- ifelse(permuted$swap,
                                  ifelse(permuted$Group == "DF1", "DF2", "DF1"),
                                  permuted$Group)

    DF1permute <- permuted[permuted$Group_perm == "DF1", , drop = FALSE]
    DF2permute <- permuted[permuted$Group_perm == "DF2", , drop = FALSE]

    # Prep and fit
    DF1run <- HMetaGroupPrep_vr(DF1permute, "DF1")
    DF2run <- HMetaGroupPrep_vr(DF2permute, "DF2")

    out1 <- try(metad_only_group(
      nR_S1 = list(as.data.frame(DF1run[[1]])),
      nR_S2 = list(as.data.frame(DF1run[[2]]))
    ), silent = TRUE)

    out2 <- try(metad_only_group(
      nR_S1 = list(as.data.frame(DF2run[[1]])),
      nR_S2 = list(as.data.frame(DF2run[[2]]))
    ), silent = TRUE)

    r.random <- matrix(NA_real_, 3, 3)

    if (!inherits(out1, "try-error") && !inherits(out2, "try-error")) {
      v1 <- summary(out1)
      v2 <- summary(out2)

      s1 <- data.frame(mean = v1$statistics[, "Mean"], name = rownames(v1$statistics))
      s2 <- data.frame(mean = v2$statistics[, "Mean"], name = rownames(v2$statistics))

      r.random[1,1] <- exp(s1$mean[s1$name == "mu_logMratio"])
      r.random[2,1] <- exp(s2$mean[s2$name == "mu_logMratio"])
      r.random[3,1] <- r.random[2,1] - r.random[1,1]

      r.random[1,2] <- abs(s1$mean[s1$name == "mu_d1"])
      r.random[2,2] <- abs(s2$mean[s2$name == "mu_d1"])
      r.random[3,2] <- r.random[2,2] - r.random[1,2]

      r.random[1,3] <- abs(s1$mean[s1$name == "mu_meta_d"])
      r.random[2,3] <- abs(s2$mean[s2$name == "mu_meta_d"])
      r.random[3,3] <- r.random[2,3] - r.random[1,3]
    }

    data.frame(
      rep          = i,
      muGroup1     = r.random[1,1],
      muGroup2     = r.random[2,1],
      mu_diff      = r.random[3,1],
      d1Group1     = r.random[1,2],
      d1Group2     = r.random[2,2],
      d1_diff      = r.random[3,2],
      meta_dGroup1 = r.random[1,3],
      meta_dGroup2 = r.random[2,3],
      meta_d_diff  = r.random[3,3],
      stringsAsFactors = FALSE
    )
  }

  # Clean up parallel backend
  doParallel::stopImplicitCluster()

  permuted_metad
}


make_results <- function(metaDat, group_label = c("CCD", "NT")[1], type_name = c("BINO", "LAT", "MONO", "ONL", "LAB")[1], simulate = 0) {
  # Map from dummy type to the short label your HMeta function expects
  type_map <- c(BINO = "BIN", LAT = "LAT", MONO = "MON", ONL = 'ONL', LAB = 'LAB')

  # Run fit for given type
  fit <- metad_only_group(
    nR_S1 = list(metaDat[[1]]),
    nR_S2 = list(metaDat[[2]]),
    sim = simulate
  )

  post <- HMeta_post_clean_meta_d(fit, paste(group_label, type_map[[type_name]]))
  df   <- post$Fit

  # Add identifiers
  df$type  <- type_name
  df$group <- group_label

  df
}

plot_permute <- function(perm_df, diff_col = 4, true_diff, p_val){

  plot_perm <- perm_df %>%
  as.data.frame() %>%
  rename(Diff = diff_col) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = true_diff + 0.15, y = 7.5, label = true_diff, fill = colpal_vr[2], colour = 'white')+
  geom_vline(xintercept = c(true_diff), color = 'black', size = 1.5) +
  geom_text(x = -1, y = 7.5, label = paste('p = ', p_val), check_overlap = T, fontface = 'bold') +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = 'cm'))

  return(plot_perm)
}

perm_effect_stats <- function(x_obs, perm) {
  # Null distribution mean & SD
  mu_null <- mean(perm)
  sd_null <- sd(perm)

  # 1) Standardized effect (null SD units)
  d_perm <- (x_obs - mu_null) / sd_null

  # 2) Probability of superiority (U)
  U <- mean(perm < x_obs)

  # 3) Two-sided permutation p-value
  p_one <- mean(perm >= abs(x_obs))

  # 4) Rank–biserial correlation
  r_rb <- 2 * U - 1

  # 5) 95% Wilson CI for U
  B <- length(perm)
  z <- qnorm(0.975)
  den <- 1 + z^2 / B
  center <- (U + z^2 / (2 * B)) / den
  half <- z * sqrt((U * (1 - U) / B) + z^2 / (4 * B^2)) / den
  U_ci <- c(lower = center - half, upper = center + half)

  list(
    d_perm = d_perm,
    U = U,
    r_rank_biserial = r_rb,
    p_one_sided = p_one,
    U_CI_95 = U_ci
  )
}

run_recovery <- function(input_data, task_label, group_label) {

  message("Running initial fit with simulation...")
  fit_initial <- metad_group(
    nR_S1 = list(input_data[[1]]),
    nR_S2 = list(input_data[[2]]),
    sim = 1
  )

  message("Extracting real parameters...")
  post <- HMeta_post_clean_meta_d(fit_initial, paste(task_label, group_label))
  df   <- post$Fit
  df$type  <- group_label
  df$group <- task_label

  real_parms <- df %>%
    filter(!str_detect(name, "counts_sim")) %>%
    rename(mean_real = mean,
           lower_real = lower,
           upper_real = upper)

  message("Extracting simulated counts...")
  sim_counts_raw <- spread_draws(fit_initial, counts_sim[subj, bin]) %>%
    group_by(subj, bin) %>%
    summarise(predicted_mean = mean(counts_sim), .groups = 'drop')

  n_subj <- length(unique(sim_counts_raw$subj))
  n_rows <- nrow(input_data[[1]])

  rec_list <- list()
  rec_list[[1]] <- matrix(NA, ncol = n_subj, nrow = n_rows)
  rec_list[[2]] <- matrix(NA, ncol = n_subj, nrow = n_rows)

  for(i in 1:n_subj) {
    subj_data <- sim_counts_raw %>%
      filter(subj == i) %>%
      arrange(bin) %>%
      pull(predicted_mean) %>%
      round()

    rec_list[[1]][,i]  <- subj_data[1:n_rows]
    rec_list[[2]][, i] <- subj_data[(n_rows + 1):(n_rows * 2)]
  }

  rec_list[[1]] <- as.data.frame(rec_list[[1]])
  rec_list[[2]] <- as.data.frame(rec_list[[2]])

  message("Running recovered fit...")
  fit_initial_rec <- metad_group(
    nR_S1 = list(rec_list[[1]]),
    nR_S2 = list(rec_list[[2]]),
    sim = 0
  )

  message("Extracting recovered fit...")
  post_rec <- HMeta_post_clean_meta_d(fit_initial_rec, paste(task_label, group_label))
  df_rec   <- post_rec$Fit
  df_rec$type  <- group_label
  df_rec$group <- task_label

  recovered_parms <- df_rec %>%
    filter(!str_detect(name, "counts_sim")) %>%
    rename(mean_rec = mean,
           lower_rec = lower,
           upper_rec = upper)

  message("Merging results...")
  final_results <- cbind(
    real_parms,
    recovered_parms %>%
      dplyr::select(mean_rec, lower_rec, upper_rec)
  )

  return(list(final_results = final_results,
              sim_data = rec_list)
         )
}


run_indiv_recovery <- function(input_data, task_label, group_label) {

  message("Running initial fit with simulation...")
  fit_initial <- metad_group(
    nR_S1 = list(input_data[[1]]),
    nR_S2 = list(input_data[[2]]),
    sim = 1
  )

  message("Extracting real parameters...")
  post <- HMeta_post_clean_meta_d(fit_initial, paste(task_label, group_label))
  real_parms <- post$Fit %>%
  filter(str_detect(name, "Mratio"),
         !str_detect(name, "log")) %>%
  rename(mean_real = mean, lower_real = lower, upper_real = upper) %>%
    mutate(type = group_label, group = task_label)

  message("Extracting simulated counts...")
  sim_counts_raw <- spread_draws(fit_initial, counts_sim[subj, bin]) %>%
    group_by(subj, bin) %>%
    summarise(predicted_mean = mean(counts_sim), .groups = 'drop')

  n_subj <- length(unique(sim_counts_raw$subj))
  n_rows <- nrow(input_data[[1]]) # Number of bins per S category

  # Storage for individual recovery results
  all_subjects_rec <- list()

  for (i in 1:n_subj) {
    message(paste("Recovering Subject", i, "of", n_subj, "..."))

    subj_full_data <- sim_counts_raw %>%
      filter(subj == i) %>%
      arrange(bin) %>%
      pull(predicted_mean) %>%
      round()

    rec_s1 <- subj_full_data[1:n_rows]
    rec_s2 <- subj_full_data[(n_rows + 1):(n_rows * 2)]

    fit_indiv_rec <- metad_indiv(nR_S1 = rec_s1, nR_S2 = rec_s2)

    post_rec <- HMeta_post_clean_meta_d(fit_indiv_rec, paste(task_label, group_label))

    all_subjects_rec[[i]] <- post_rec$Fit %>%
      filter(!str_detect(name, "counts_sim")) %>%
      dplyr::select(mean_rec = mean, lower_rec = lower, upper_rec = upper) %>%
      mutate(subj_id = i)%>%
      slice(2)
  }

  rec_df_final  <- bind_rows(all_subjects_rec)

  final_results <- cbind(real_parms, rec_df_final)

  return(final_results)
}


make_vr_perm_plot <- function(perm_df, diff_val) {
  perm_df %>% as.data.frame() %>% rename(Diff = 4) %>%
    ggplot(aes(Diff)) + geom_histogram(binwidth = 0.01) +
    labs(x = expression(paste("Randomised ",mu['ratio'])), y = 'Sample Count') +
    geom_label(x = diff_val + 0.15, y = 7.5, label = round(diff_val,2),
               fill = colpal_vr[2], colour = 'white') +
    geom_vline(xintercept = diff_val, color = colpal_vr[2], size = 1.5) +
    theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=24))
}

# Decision-phase RTs (confidence not reported; conf == 0)
run_rt_lmer <- function(data, env_filter, conf_filter) {
  lmerTest::lmer(log(RT) ~ CCD * kmed + (1|ID),
                 data = data %>% filter(type=='main', conf==conf_filter,
                                        kmed!='kmed', environment==env_filter)) %>% summary()
}

# Coherence trajectories for Experiment 3 (each presentation type)
make_vr_cal_plot <- function(ptype, legend_name) {
    ggplot(check_coh_RDK %>%
             filter(type == 'Calibration', prestype == ptype, vr == 1) %>%
             group_by(ID) %>% mutate(Trial = 1:n())) +
      stat_summary(aes(Trial, coherence, fill = group), geom='ribbon', alpha=0.1, colour=NA) +
      stat_summary(aes(Trial, coherence, colour = group), geom='line') +
      scale_color_manual(values = colpal_vr, name = legend_name) +
      scale_fill_manual(values  = colpal_vr, name = legend_name) +
      labs(x = 'Trial', y = 'Coherence') +
      coord_cartesian(ylim = c(0, 0.5), xlim = c(0, ifelse(ptype=='Binocular', 24, 48))) +
      scale_x_continuous(expand = c(0,0)) +
      theme_minimal(base_size = 20) +
      theme(legend.position = 'top', panel.grid.minor = element_blank())
}

# Experiment 3: per-presentation-type beeswarm panels
make_exp3_plot <- function(data, y_var, y_lab, label_y) {
  lapply(c('Bino.', 'Lat.', 'Mono.'), function(p) {
    tidyplot(data %>% filter(Pres == p), x = group, y = !!sym(y_var), colour = kmed) %>%
      add_data_points_beeswarm(alpha=0.1, cex=1) %>% add_sd_errorbar() %>% add_mean_dot(size=2) %>%
      add_title(p) %>%
      adjust_size(width=75, height=50) %>% adjust_x_axis_title(title='') %>%
      adjust_y_axis_title(title=y_lab) %>% adjust_font(fontsize=15) %>%
      adjust_legend_title(title='Coh.') %>% adjust_legend_position(position='top') +
      stat_compare_means(paired=TRUE, hide.ns=TRUE, size=5, show.legend=FALSE,
                         label='p.signif', label.y=label_y) +
      coord_cartesian(clip="off")
  }) %>% Reduce(`|`, .)
}


# Helper to reformat a matrix into long format for per-participant correlation
fit_to_count <- function(mat, type_label, exp_label, dec_label) {
  as.data.frame(mat) %>%
    mutate(Bin = 1:n()) %>%
    pivot_longer(cols = starts_with("V"), names_to = "Participant", values_to = "Count") %>%
    mutate(Type = type_label, Exp = exp_label, Dec = dec_label)
}

make_vr_hist <- function(groups_to_plot) {
  mcmcGroup_vr %>%
    filter(Parameter == "mu_logMratio", group %in% groups_to_plot) %>%
    ggplot(aes(exp(value), fill = group)) +
    geom_histogram(binwidth = 0.03, alpha = 0.9, position = 'dodge') +
    scale_fill_manual(values = colpal_vr[1:2]) +
    labs(y = "Sample Count", x = expression(paste(mu["ratio"]))) +
    hmeta_theme +
    theme(legend.position = 'none')
}

# Helper to run model + effect sizes in one call
run_lmer <- function(formula, data) {
  m <- lmerTest::lmer(formula, data = data)
  e <- effectsize::effectsize(m, method = "posthoc")
  return(list(summary = summary(m), effect_size = e))
}

# Summary table of group-level mu and sigma for mu_logMratio (exp-transformed
# to give the M-ratio on its natural scale)
extract_stat <- function(res, param) {
  data.frame(
    mean  = exp(res[[3]]$mean[ res[[3]]$name == param]),
    lower = exp(res[[3]]$lower[res[[3]]$name == param]),
    upper = exp(res[[3]]$upper[res[[3]]$name == param])
  )
}

# Reshape ccd_b / ccd_l / ccd_m etc. into a long format aligned by participant
make_block <- function(group, type, v) {
  data.frame(group = group, type = type, value = v, stringsAsFactors = FALSE)
}
