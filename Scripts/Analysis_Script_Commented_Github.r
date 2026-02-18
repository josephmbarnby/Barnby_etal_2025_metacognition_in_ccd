
# Barnby et al. (2026): Corpus Callosum Dysgenesis impairs metacognition:
# evidence from multi-modality and multi-cohort replications
#
# Analysis script covering three experiments:
#   Experiment 1: Online RDK (CCD vs NT, computer-based)
#   Experiment 2: In-lab RDK (CCD vs NT, MRI environment)
#   Experiment 3: VR RDK (CCD vs NT, binocular / monocular / lateralized)
#
# Key outcome measures: perceptual accuracy, confidence, metacognitive
# efficiency (hierarchical meta-d' / d'; HMeta-D method, Maniscalco & Lau 2012)

# ── 1. SETUP ------------------------------------------------------------

library(ggpubr);    library(tidyverse);  library(patchwork)
library(doParallel);library(rstanarm);   library(brms)
library(magrittr);  library(reshape2);   library(rjags)
library(coda);      library(lattice);    library(lme4)
library(broom);     library(ggmcmc);     library(foreach)
library(sjPlot);    library(webshot);    library(tidyplots); library(tidybayes)

source('RDKUtilityFuncs.R')   # custom helper functions (coherence prep,
                               # HMeta formatting, permutation wrappers, etc.)

# Colour palettes used throughout all figures
colpal    <- c('#931F1D', '#FE938C', '#475B5A', '#81D6E3')  # Exp 1 & 2 (4 groups)
colpal_vr <- c('#B80C09', '#0E1C36')                        # Exp 3 (CCD / NT)

# ── 2. LOAD DATA ------------------------------------------------------------

# Experiments 1 (online) and 2 (MRI/lab): conventional RDK task
Experiment_1_CCD <- read.csv('Data/Experiment1_CCD.csv')
Experiment_1_NT  <- read.csv('Data/Experiment1_NT.csv')
Experiment_2_CCD <- read.csv('Data/Experiment2_CCD.csv')
Experiment_2_NT  <- read.csv('Data/Experiment2_NT.csv')

# Experiment 3: VR task with binocular / monocular / lateralized presentations
data_vr <- read.csv('Data/Experiment3.csv')

# ── 3. DATA QUALITY: CALIBRATION & EXCLUSIONS ------------------------------------------------------------
# Participants were excluded if their coherence threshold during calibration
# exceeded 0.5, if confidence ratings showed no variance, or if self-reported
# error counts were implausibly high. Trials on which participants reported a
# button-press mistake are also removed.

## --- 3a. Visual diagnostics ----
# Inspect coherence trajectories across calibration trials to flag outliers
(Experiment_1_NT  %>% check_coherence('NT (Online)', legend = F) +
 Experiment_1_CCD %>% check_coherence('CCD (Online)', legend = F)) /
(Experiment_2_NT  %>% check_coherence('NT (MRI)',    legend = F) +
 Experiment_2_CCD %>% check_coherence('CCD (MRI)',   legend = F))

# Spot-check specific outlier participants identified during calibration review
Experiment_2_CCD %>%
  filter(ID %in% c('1001', '1027')) %>%
  check_coherence('CCD (MRI)', legend = F)

# Inspect confidence-rating distributions per participant (ridge plots)
# — applied to each of the four groups; code is identical so only one shown
for (dat in list(Experiment_1_CCD, Experiment_2_CCD, Experiment_1_NT, Experiment_2_NT)) {
  P_dat <- dat %>%
    filter(conf != 0, type == 'main', mistake %in% c(0, NA, FALSE)) %>%
    ggplot(aes(conf, factor(ID), fill = ID, group = ID)) +
    ggridges::geom_density_ridges(scale = 1, jittered_points = TRUE,
                                  position = "raincloud", alpha = 0.7, scale = 0.2) +
    labs(x = 'Confidence', y = 'ID') +
    theme_bw() +
    theme(text = element_text(size = 18), legend.position = 'none')
  print(P_dat)
}

# Count self-reported button-press errors per participant (used to set
# the ">20 mistakes" exclusion threshold applied below)
for (dat in list(Experiment_1_CCD, Experiment_2_CCD, Experiment_1_NT, Experiment_2_NT)) {
  P_dat <- dat %>%
    filter(conf != 0, type == 'main') %>%
    group_by(ID) %>%
    mutate(mistake_sum = ifelse(!is.na(mistake) & mistake == TRUE, 1, 0),
           mistake_sum = sum(mistake_sum)) %>%
    dplyr::select(mistake_sum, ID) %>%
    unique() %>%
    ggplot(aes(factor(ID), mistake_sum, fill = ID)) +
    geom_col() +
    coord_cartesian(ylim = c(1, 25)) +
    labs(x = 'ID', y = 'Mistake') +
    theme_bw() +
    theme(text = element_text(size = 18), legend.position = 'none')
  print(P_dat)
}

## --- 3b. Apply exclusions ----
# Reasons documented per participant:
#   2004  : coherence threshold > 0.5 (Exp 1 CCD)
#   5f20… : coherence threshold > 0.5 (Exp 1 NT)
#   1001  : coherence threshold > 0.5 (Exp 2 CCD)
#   1006  : zero variance in confidence + >20 mistakes (Exp 2 CCD)
#   1027  : >20 mistakes (Exp 2 CCD)
#   3002  : zero variance in confidence (Exp 2 NT)
#   1004  : self-reported mistakes not marked in data (Exp 2 CCD)
#   2008  : uncertain phenotype (Exp 1 CCD); retained in sensitivity set

Experiment_1_CCDcheck       <- Experiment_1_CCD %>% filter(ID != '2004')
Experiment_1_CCDcheck_w_exc <- Experiment_1_CCD %>% filter(ID != '2004', ID != '2008')
Experiment_1_NTcheck        <- Experiment_1_NT  %>% filter(ID != '5f20a850b6d5ed3d186bc3c3')
Experiment_2_CCDcheck       <- Experiment_2_CCD %>% filter(!ID %in% c('1001', '1006', '1027'))
Experiment_2_NTcheck        <- Experiment_2_NT  %>% filter(ID != '3002')

# Remove mistake-flagged trials across all groups
Experiment_1_CCDcheck <- Experiment_1_CCDcheck %>% filter(mistake == FALSE)
Experiment_1_NTcheck  <- Experiment_1_NTcheck  %>% filter(mistake == FALSE)
Experiment_2_CCDcheck <- Experiment_2_CCDcheck %>% filter(mistake == FALSE)
Experiment_2_NTcheck  <- Experiment_2_NTcheck  %>% filter(mistake == FALSE)

# Final sample sizes (Experiments 1 & 2)
cat("Exp1 CCD (main):", length(unique(Experiment_1_CCDcheck$ID)),       "\n")
cat("Exp1 CCD (w/ 2008 excluded):", length(unique(Experiment_1_CCDcheck_w_exc$ID)), "\n")
cat("Exp1 NT:",         length(unique(Experiment_1_NTcheck$ID)),         "\n")
cat("Exp2 CCD:",        length(unique(Experiment_2_CCDcheck$ID)),         "\n")
cat("Exp2 NT:",         length(unique(Experiment_2_NTcheck$ID)),          "\n")

checkBothRDK <- rbind(
  Experiment_1_CCDcheck  %>% mutate(correct = as.numeric(correct)),
  Experiment_1_NTcheck   %>% mutate(correct = as.numeric(correct)),
  Experiment_2_CCDcheck  %>% mutate(correct = as.numeric(correct)),
  Experiment_2_NTcheck   %>% mutate(correct = as.numeric(correct))
)

# Sensitivity dataset (Exp 1 CCD with participant 2008 also excluded)
checkBothRDK_w_exc <- rbind(
  Experiment_1_CCDcheck_w_exc %>% mutate(correct = as.numeric(correct)),
  Experiment_1_NTcheck        %>% mutate(correct = as.numeric(correct)),
  Experiment_2_CCDcheck       %>% mutate(correct = as.numeric(correct)),
  Experiment_2_NTcheck        %>% mutate(correct = as.numeric(correct))
)

# Binary flag for the high-coherence (kmed × 2) condition
checkBothRDK$highCoh       <- ifelse(checkBothRDK$kmed       == "kmed x 2", 1, 0)
checkBothRDK_w_exc$highCoh <- ifelse(checkBothRDK_w_exc$kmed == "kmed x 2", 1, 0)

# ── 4. ENVIRONMENT CHECK & DEMOGRAPHICS -------------------------
# We test whether the testing platform (computer / MRI / VR) explains variance
# in accuracy or confidence. If not, experiments can be analysed jointly.

# Accuracy by environment
accuracy_data_vr <- data_vr %>%
  rename(ID = PID, environment = Environment) %>%
  filter(TrialType == 'Main', Presentation == 'Binocular') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Accuracy = mean(Correct), .groups = 'drop')
accuracy_data_vr$ID <- as.character(accuracy_data_vr$ID)

accuracy_data_non_vr <- checkBothRDK %>%
  filter(type == 'main') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Accuracy = mean(correct), .groups = 'drop')

accuracy_data <- rbind(accuracy_data_vr, accuracy_data_non_vr)

# Linear model: environment and CCD group as predictors of accuracy
model_acc <- lm(Accuracy ~ environment + CCD, data = accuracy_data)
summary(aov(Accuracy ~ environment + CCD, data = accuracy_data))
plot_model(model_acc, title = "Effect of environment on participant accuracy",
           show.values = TRUE, show.p = TRUE)
tab_model(model_acc, show.re.var = TRUE,
          dv.labels = "Effect of environment on participant accuracy",
          file = "Stats/participant_accuracy_lmer.html")

# Confidence by environment (scale confidence to z-scores before aggregating)
data_vr$Confidence <- scale(data_vr$ReportedVeryConfident)
confidence_data_vr <- data_vr %>%
  filter(TrialType == 'Main') %>%
  group_by(PID, Environment, CCD) %>%
  summarise(Confidence = mean(Confidence), .groups = 'drop') %>%
  mutate(PID = as.character(PID)) %>%
  rename(ID = PID, environment = Environment)

checkBothRDKConfidenceTrialsOnly <- checkBothRDK %>%
  filter(conf != 0) %>%
  mutate(Confidence = scale(conf))

confidence_data_non_vr <- checkBothRDKConfidenceTrialsOnly %>%
  filter(type == 'main') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Confidence = mean(Confidence), .groups = 'drop')

confidence_data <- rbind(confidence_data_vr, confidence_data_non_vr)
summary(aov(Confidence ~ environment + CCD, data = confidence_data))

# ── 5. BEHAVIOURAL RESULTS: EXPERIMENTS 1 & 2 -----------------------------------
# Primary outcome: does accuracy / confidence scale with coherence level,
# and does this scaling differ between CCD and NT?

lmer_data     <- checkBothRDK       %>% filter(type == 'main', kmed != 'kmed', mistake %in% c(0,NA,FALSE), conf != 0)
lmer_data_exc <- checkBothRDK_w_exc %>% filter(type == 'main', kmed != 'kmed', mistake %in% c(0,NA,FALSE), conf != 0)
data_vr_tests <- data_vr %>%
  rename(ID=PID, type=TrialType, coherence=ActiveCoherence, group=Group) %>%
  filter(type=='Main') %>%
  group_by(ID, PresentationType) %>%
  mutate(
         conf = scale(ReportedVeryConfident),
         #conf = 1/(1+exp(-scale(conf))),
         kmed = ifelse(HighCoherenceTrial==1, 'x2','x0.5'),
         kmed = factor(kmed, levels = c('x0.5', 'x2'), ordered = T)
         ) %>%
  group_by(ID, group)

# --- 5a. Reaction-time density plots (sanity check) ----
CCDcheck1 <- ggplot(
  checkBothRDK %>%
    filter(type == 'main', kmed != 'kmed', RT < 5000, mistake %in% c(0,NA,FALSE)) %>%
    mutate(correct = ifelse(correct == 1, 'Correct: Yes', 'Correct: No')),
  aes(RT, fill = group)) +
  geom_density(show.legend = TRUE, alpha = 0.7) +
  scale_fill_manual(values = colpal) +
  facet_wrap(kmed ~ correct) +
  labs(x = 'Reaction Time (ms)') +
  theme_minimal() +
  theme(axis.title.y = element_blank())

# Non-parametric tests on RT medians (group × environment)
kruskal.test(RT ~ group, data = rbind(Experiment_1_CCDcheck, Experiment_1_NTcheck))
kruskal.test(RT ~ group, data = rbind(Experiment_2_CCDcheck, Experiment_2_NTcheck))
kruskal.test(RT ~ group, data = rbind(Experiment_2_CCDcheck, Experiment_1_CCDcheck))
kruskal.test(RT ~ group, data = rbind(Experiment_1_NTcheck,  Experiment_1_NTcheck))

# --- 5b. Per-participant accuracy and confidence by coherence level --------
# Compute individual means across the two coherence conditions

## EXPERIMENT 1 & 2 -------

av_cor <- checkBothRDK %>%
  filter(type == 'main', kmed != 'kmed', mistake %in% c(0,NA,FALSE)) %>%
  group_by(group, ID, kmed) %>%
  mutate(av_cor = mean(corAdjusted),
         kmed   = ifelse(kmed == 'kmed x 2', 'x2', 'x0.5')) %>%
  group_by(group, kmed) %>%
  dplyr::select(ID, av_cor, kmed) %>%
  distinct()

av_co <- checkBothRDK %>%
  filter(type == 'main', kmed != 'kmed', conf != 0, mistake %in% c(0,NA,FALSE)) %>%
  group_by(group, ID, kmed) %>%
  mutate(av_co = mean(conf),
         kmed  = ifelse(kmed == 'kmed x 2', 'x2', 'x0.5')) %>%
  group_by(group, kmed) %>%
  dplyr::select(ID, av_co, kmed) %>%
  distinct()

# Shared aesthetics for accuracy and confidence plots
dot_layer <- list(
  geom_jitter(shape = 21, width = 0.1, alpha = 0.1, colour = 'black'),
  stat_summary(),
  stat_summary(geom = 'line', show.legend = FALSE),
  scale_color_manual(values = colpal),
  scale_fill_manual(values = colpal),
  theme_minimal(),
  theme(panel.grid = element_blank(), text = element_text(size = 18))
)

CCDcheck3 <- ggplot(av_cor, aes(kmed, av_cor, color = group, group = group,
                                 fill = group)) +
  dot_layer +
  labs(x = 'Coherence', y = 'Correct') +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(legend.position = 'none')

CCDcheck4 <- ggplot(av_co, aes(kmed, av_co, color = group, group = group,
                                fill = group)) +
  dot_layer +
  labs(x = 'Coherence', y = 'Confidence') +
  coord_cartesian(ylim = c(50, 100)) +
  theme(legend.position = 'none')

(CCDcheck3 | CCDcheck4)

## EXPERIMENT 3 --------

# Per-participant accuracy and confidence by coherence × presentation type
av_cor_vr <- data_vr_tests %>%
  group_by(group, ID, kmed, PresentationType) %>%
  summarise(av_cor = mean(Correct), .groups='drop')

av_co_vr <- data_vr_tests %>%
  group_by(group, ID, kmed, PresentationType) %>%
  summarise(av_co = mean(conf), .groups='drop')

# Accuracy and confidence plots (faceted by presentation type)
vr_cor_plot <- ggplot(av_cor_vr, aes(kmed, av_cor, color = group, group = group, fill = group)) +
  geom_jitter(shape = 21, width = 0.1, alpha = 0.1, colour = 'black') +
  stat_summary() + stat_summary(geom = 'line', show.legend = FALSE) +
  labs(x = 'Coherence', y = 'Correct') +
  facet_wrap(~PresentationType) +
  scale_color_manual(values = colpal_vr) + scale_fill_manual(values = colpal_vr) +
  coord_cartesian(ylim = c(0.5, 1)) + theme_minimal() +
  theme(legend.position = 'none', axis.title.x = element_blank())

vr_co_plot <- ggplot(av_co_vr, aes(kmed, av_co, color = group, group = group, fill = group)) +
  geom_jitter(shape = 21, width = 0.1, alpha = 0.1, colour = 'black', show.legend = FALSE) +
  stat_summary() + stat_summary(geom = 'line', show.legend = FALSE) +
  labs(x = 'Coherence', y = 'Confidence') +
  facet_wrap(~PresentationType) +
  scale_color_manual(values = colpal_vr, name = 'Group') +
  scale_fill_manual(values = colpal_vr) +
  theme_minimal() +
  theme(legend.position = 'none', strip.text.x = element_blank())

(vr_cor_plot / vr_co_plot) & theme(text = element_text(size = 18), panel.grid = element_blank())

# --- 5c. Paired t-tests: coherence effect within each group ----------

for (grp in c('CCD (Online)', 'NT (Online)', 'CCD (MRI)', 'NT (MRI)')) {
  cat("\n--- Accuracy:", grp, "---\n")
  print(t.test(av_cor ~ kmed, data = av_cor %>% filter(group == grp), paired = TRUE))
  cat("\n--- Confidence:", grp, "---\n")
  print(t.test(av_co  ~ kmed, data = av_co  %>% filter(group == grp), paired = TRUE))
}

# Group × coherence ANOVAs (Exp 1 and Exp 2 separately)
aov(av_cor ~ group*kmed, data = av_cor %>%
      filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
aov(av_cor ~ group*kmed, data = av_cor %>%
      filter(group %in% c('NT (MRI)',    'CCD (MRI)')))    %>% summary()

aov(av_co ~ group*kmed, data = av_co %>%
      filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
aov(av_co ~ group*kmed, data = av_co %>%
      filter(group %in% c('NT (MRI)',    'CCD (MRI)')))    %>% summary()

# --- 5d. Linear mixed-effects models (trial-level; primary inference) --------
# kmed × group interaction on accuracy and confidence, with random intercept
# per participant. Sensitivity models include/exclude participant 2008.

### Experiment 1 ------

run_lmer(correct ~ kmed + (1|ID), lmer_data     %>% filter(group == 'CCD (Online)'))
run_lmer(correct ~ kmed + (1|ID), lmer_data     %>% filter(group == 'NT (Online)'))
run_lmer(correct ~ kmed * group + (1|ID), lmer_data     %>% filter(group %in% c('NT (Online)', 'CCD (Online)')))

run_lmer(conf ~ kmed + (1|ID), lmer_data     %>% filter(group == 'CCD (Online)'))
run_lmer(conf ~ kmed + (1|ID), lmer_data     %>% filter(group == 'NT (Online)'))
run_lmer(conf ~ kmed * group + (1|ID), lmer_data     %>% filter(group %in% c('NT (Online)', 'CCD (Online)')))

### Experiment 2 --------

run_lmer(correct ~ kmed + (1|ID), lmer_data %>% filter(group == 'CCD (MRI)'))
run_lmer(correct ~ kmed + (1|ID), lmer_data %>% filter(group == 'NT (MRI)'))
run_lmer(correct ~ kmed * group + (1|ID), lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)')))

run_lmer(conf ~ kmed + (1|ID), lmer_data %>% filter(group == 'CCD (MRI)'))
run_lmer(conf ~ kmed + (1|ID), lmer_data %>% filter(group == 'NT (MRI)'))
run_lmer(conf ~ kmed * group + (1|ID), lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)')))

### Experiment 3 ---------

run_lmer(Correct  ~ kmed * group + (1|ID), data_vr_tests)
run_lmer(conf  ~ kmed * group + (1|ID), data_vr_tests)

#### BINO ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular'))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular'))

run_lmer(Correct  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Binocular'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular'))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular'))

run_lmer(conf  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Binocular'))

#### LAT ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized'))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized'))

run_lmer(Correct  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Lateralized'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized'))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized'))

run_lmer(conf  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Lateralized'))

#### MONO ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular'))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular'))

run_lmer(Correct  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Monocular'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular'))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular'))

run_lmer(conf  ~ kmed * group + (1|ID), data_vr_tests %>% filter(PresentationType=='Monocular'))

# --- 5e. Confidence Bin Spread --------

LAB_CCD_metaDat    <- HMetaGroupPrep(checkBothRDK, 'CCD (MRI)')
LAB_CCD_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'CCD (MRI)')
LAB_CCD_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'CCD (MRI)')

# Create data set for ONL CCD participants
ONL_CCD_metaDat    <- HMetaGroupPrep(checkBothRDK, 'CCD (Online)')
ONL_CCD_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'CCD (Online)')
ONL_CCD_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'CCD (Online)')

# Create data set for ONL NT participants
ONL_NT_metaDat    <- HMetaGroupPrep(checkBothRDK, 'NT (Online)')
ONL_NT_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'NT (Online)')
ONL_NT_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'NT (Online)')

# Create data set for LAB CCD participants
LAB_NT_metaDat    <- HMetaGroupPrep(checkBothRDK,  'NT (MRI)')
LAB_NT_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'NT (MRI)')
LAB_NT_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'NT (MRI)')

# Confidence spread (normalised response count distributions by coherence)
conf_spread <- data.frame(
  conf = c(6:1, -1:-6),
  cond = c(rep('lab', 72), rep('online', 72)),
  type = rep(c(rep('x2',24), rep('Both',24), rep('x0.5',24)), 2),
  s1   = c(rowSums(LAB_CCD_metaDat_hc$S1mat), rowSums(LAB_NT_metaDat_hc$S1mat),
           rowSums(LAB_CCD_metaDat$S1mat),    rowSums(LAB_NT_metaDat$S1mat),
           rowSums(LAB_CCD_metaDat_lc$S1mat), rowSums(LAB_NT_metaDat_lc$S1mat),
           rowSums(ONL_CCD_metaDat_hc$S1mat), rowSums(ONL_NT_metaDat_hc$S1mat),
           rowSums(ONL_CCD_metaDat$S1mat),    rowSums(ONL_NT_metaDat$S1mat),
           rowSums(ONL_CCD_metaDat_lc$S1mat), rowSums(ONL_NT_metaDat_lc$S1mat)),
  s2   = c(rev(rowSums(LAB_CCD_metaDat_hc$S2mat)), rev(rowSums(LAB_NT_metaDat_hc$S2mat)),
           rev(rowSums(LAB_CCD_metaDat$S2mat)),     rev(rowSums(LAB_NT_metaDat$S2mat)),
           rev(rowSums(LAB_NT_metaDat_lc$S2mat)),   rev(rowSums(LAB_NT_metaDat_lc$S2mat)),
           rev(rowSums(ONL_CCD_metaDat_hc$S2mat)),  rev(rowSums(ONL_NT_metaDat_hc$S2mat)),
           rev(rowSums(ONL_CCD_metaDat$S2mat)),      rev(rowSums(ONL_NT_metaDat$S2mat)),
           rev(rowSums(ONL_NT_metaDat_lc$S2mat)),   rev(rowSums(ONL_NT_metaDat_lc$S2mat)))
) %>%
  pivot_longer(s1:s2, names_to = 'Dec', values_to = 'Count') %>%
  mutate(group = c(rep('CCD',144), rep('NT',144))) %>%
  group_by(type, Dec, cond, group) %>%
  mutate(type  = factor(type, levels = c('x0.5','Both','x2')),
         Count = Count / sum(Count)) %>%
  ggplot(aes(conf, Count, fill = group)) +
  geom_col(position = 'dodge') +
  labs(x = 'Confidence Bin (1=low, 6=high)', y = 'Normalised Count') +
  scale_x_continuous(breaks = c(6,1,-1,-6)) +
  facet_wrap(~type) +
  theme_bw(base_size = 24) +
  theme(legend.position = c(0.10,0.8), legend.background = element_blank(),
        legend.title = element_blank(), panel.grid = element_blank())

conf_spread

# ── 6. HIERARCHICAL META-D' FITTING ---------------
# Metacognitive efficiency (meta-d'/d') estimated using the hierarchical Bayesian
# method (HMeta-D; Fleming 2017). Response matrices are built with HMetaGroupPrep()
# then passed to metad_group(), which is slow and should only be run once.

# --- 6a. Prepare response count matrices ----
LAB_CCD_metaDat       <- HMetaGroupPrep(checkBothRDK, 'CCD (MRI)')
ONL_CCD_metaDat       <- HMetaGroupPrep(checkBothRDK, 'CCD (Online)')
ONL_CCD_metaDat_exc   <- HMetaGroupPrep(checkBothRDK_w_exc, 'CCD (Online)')

ONL_NT_metaDat        <- HMetaGroupPrep(checkBothRDK, 'NT (Online)')
LAB_NT_metaDat        <- HMetaGroupPrep(checkBothRDK, 'NT (MRI)')

# VR: three presentation types × with/without sensitivity exclusions
VR_CCD_metaDat        <- HMetaGroupPrep_vr(data_vr %>% filter(BinocularTrial   == 1), 'CCD')
VR_CCD_metaDat_L      <- HMetaGroupPrep_vr(data_vr %>% filter(LateralizedTrial == 1), 'CCD')
VR_CCD_metaDat_M      <- HMetaGroupPrep_vr(data_vr %>% filter(MonocularTrial   == 1), 'CCD')
VR_CCD_metaDat_exc    <- HMetaGroupPrep_vr(data_vr %>% filter(BinocularTrial   == 1, !PID %in% exc_ids), 'CCD')
VR_CCD_metaDat_L_exc  <- HMetaGroupPrep_vr(data_vr %>% filter(LateralizedTrial == 1, !PID %in% exc_ids), 'CCD')
VR_CCD_metaDat_M_exc  <- HMetaGroupPrep_vr(data_vr %>% filter(MonocularTrial   == 1, !PID %in% exc_ids), 'CCD')

VR_NT_metaDat         <- HMetaGroupPrep_vr(data_vr %>% filter(BinocularTrial   == 1), 'NT')
VR_NT_metaDat_L       <- HMetaGroupPrep_vr(data_vr %>% filter(LateralizedTrial == 1), 'NT')
VR_NT_metaDat_M       <- HMetaGroupPrep_vr(data_vr %>% filter(MonocularTrial   == 1), 'NT')

# --- 6b. Fit HMeta-D (run once; results saved to disk) --------
## NOTE: outputs are stochastic. Expect minor numeric deviation from published work.

do_mratio <- 0
if(do_mratio == 1){

  setwd('HMeta-d-master/R/')
  source('Function_metad_group.R')

  outputOnline_NT  <- metad_group(nR_S1 = list(ONL_NT_metaDat[[1]]),  nR_S2 = list(ONL_NT_metaDat[[2]]));
  outputLab_NT     <- metad_group(nR_S1 = list(LAB_NT_metaDat[[1]]),  nR_S2 = list(LAB_NT_metaDat[[2]]));

  outputOnline_CCD  <- metad_group(nR_S1 = list(ONL_CCD_metaDat[[1]]), nR_S2 = list(ONL_CCD_metaDat[[2]]));
  outputLab_CCD     <- metad_group(nR_S1 = list(LAB_CCD_metaDat[[1]]), nR_S2 = list(LAB_CCD_metaDat[[2]]));

  outputVR_NT      <- metad_group(nR_S1 = list(VR_NT_metaDat[[1]]),   nR_S2 = list(VR_NT_metaDat[[2]]));
  outputVR_NT_L    <- metad_group(nR_S1 = list(VR_NT_metaDat_L[[1]]), nR_S2 = list(VR_NT_metaDat_L[[2]]));
  outputVR_NT_M    <- metad_group(nR_S1 = list(VR_NT_metaDat_M[[1]]), nR_S2 = list(VR_NT_metaDat_M[[2]]));

  outputVR_CCD     <- metad_group(nR_S1 = list(VR_CCD_metaDat[[1]]),  nR_S2 = list(VR_CCD_metaDat[[2]]));
  outputVR_CCD_L   <- metad_group(nR_S1 = list(VR_CCD_metaDat_L[[1]]),nR_S2 = list(VR_CCD_metaDat_L[[2]]));
  outputVR_CCD_M   <- metad_group(nR_S1 = list(VR_CCD_metaDat_M[[1]]),nR_S2 = list(VR_CCD_metaDat_M[[2]]));

  outputVR_CCD_exc      <- metad_group(nR_S1 = list(VR_CCD_metaDat_exc[[1]]),   nR_S2 = list(VR_CCD_metaDat_exc[[2]]));
  outputVR_CCD_L_exc    <- metad_group(nR_S1 = list(VR_CCD_metaDat_L_exc[[1]]), nR_S2 = list(VR_CCD_metaDat_L_exc[[2]]));
  outputVR_CCD_M_exc    <- metad_group(nR_S1 = list(VR_CCD_metaDat_M_exc[[1]]), nR_S2 = list(VR_CCD_metaDat_M_exc[[2]]));

  # Save the output
  saveRDS(outputLab_CCD,    'HMetaDFit/Mratio/LAB_CCD.rdata')
  saveRDS(outputOnline_CCD,    'HMetaDFit/Mratio/ONL_CCD.rdata')
  saveRDS(outputOnline_CCD_exc,    'HMetaDFit/Mratio/ONL_CCD_exc.rdata')
  saveRDS(outputVR_CCD,    'HMetaDFit/Mratio/VR_CCD.rdata')
  saveRDS(outputVR_CCD_L,  'HMetaDFit/Mratio/VR_CCD_L.rdata')
  saveRDS(outputVR_CCD_M,  'HMetaDFit/Mratio/VR_CCD_M.rdata')
  saveRDS(outputVR_CCD_exc,    'HMetaDFit/Mratio/VR_CCD_exc.rdata')
  saveRDS(outputVR_CCD_L_exc,  'HMetaDFit/Mratio/VR_CCD_L_exc.rdata')
  saveRDS(outputVR_CCD_M_exc,  'HMetaDFit/Mratio/VR_CCD_M_exc.rdata')
  saveRDS(outputLab_NT,    'HMetaDFit/Mratio/LAB_NT.rdata')
  saveRDS(outputOnline_NT, 'HMetaDFit/Mratio/ONL_NT.rdata')
  saveRDS(outputVR_NT,     'HMetaDFit/Mratio/VR_NT.rdata')
  saveRDS(outputVR_NT_L,   'HMetaDFit/Mratio/VR_NT_L.rdata')
  saveRDS(outputVR_NT_M,   'HMetaDFit/Mratio/VR_NT_M.rdata')

  setwd('~')

}

# --- 6c. Load pre-fitted models ----
outputLAB_CCD     <- readRDS('HMetaDFit/Mratio/LAB_CCD.rdata')
outputONL_CCD     <- readRDS('HMetaDFit/Mratio/ONL_CCD.rdata')
outputVR_CCD      <- readRDS('HMetaDFit/Mratio/VR_CCD.rdata')
outputVR_CCD_L    <- readRDS('HMetaDFit/Mratio/VR_CCD_L.rdata')
outputVR_CCD_M    <- readRDS('HMetaDFit/Mratio/VR_CCD_M.rdata')
outputLAB_NT      <- readRDS('HMetaDFit/Mratio/LAB_NT.rdata')
outputONL_NT      <- readRDS('HMetaDFit/Mratio/ONL_NT.rdata')
outputVR_NT       <- readRDS('HMetaDFit/Mratio/VR_NT.rdata')
outputVR_NT_L     <- readRDS('HMetaDFit/Mratio/VR_NT_L.rdata')
outputVR_NT_M     <- readRDS('HMetaDFit/Mratio/VR_NT_M.rdata')

# --- 6d. Extract cleaned MCMC samples and posterior summaries ----
Results_LAB_CCD     <- HMeta_post_clean(outputLAB_CCD,     'CCD (LAB)', F)
Results_ONL_CCD     <- HMeta_post_clean(outputONL_CCD,     'CCD (ONL)', F)
Results_ONL_CCD_exc <- HMeta_post_clean(outputONL_CCD_exc, 'CCD (ONL)', F)
Results_VR_CCD      <- HMeta_post_clean(outputVR_CCD,      'CCD (VR)',  F)
Results_VR_CCD_L    <- HMeta_post_clean(outputVR_CCD_L,    'CCD (VR L)',F)
Results_VR_CCD_M    <- HMeta_post_clean(outputVR_CCD_M,    'CCD (VR M)',F)
Results_ONL_NT      <- HMeta_post_clean(outputONL_NT,      'NT (ONL)',  F)
Results_LAB_NT      <- HMeta_post_clean(outputLAB_NT,      'NT (LAB)',  F)
Results_VR_NT       <- HMeta_post_clean(outputVR_NT,       'NT (VR)',   F)
Results_VR_NT_L     <- HMeta_post_clean(outputVR_NT_L,     'NT (VR L)', F)
Results_VR_NT_M     <- HMeta_post_clean(outputVR_NT_M,     'NT (VR M)', F)

# Stack MCMC draws for Experiments 1 & 2
mcmcGroup <- rbind(Results_LAB_CCD[[1]], Results_ONL_CCD[[1]],
                   Results_ONL_NT[[1]],  Results_LAB_NT[[1]])

# Stack MCMC draws for Experiment 3 (all three presentation types)
mcmcGroup_vr <- rbind(Results_VR_CCD[[1]],   Results_VR_NT[[1]],
                       Results_VR_CCD_L[[1]], Results_VR_NT_L[[1]],
                       Results_VR_CCD_M[[1]], Results_VR_NT_M[[1]])

statGroup <- data.frame(
  CCD_EXP1 = c(unlist(extract_stat(Results_ONL_CCD, "mu_logMratio")),
              unlist(extract_stat(Results_ONL_CCD, "sigma_logMratio"))),
  NT_EXP1  = c(unlist(extract_stat(Results_ONL_NT,  "mu_logMratio")),
              unlist(extract_stat(Results_ONL_NT,  "sigma_logMratio"))),
  CCD_EXP2 = c(unlist(extract_stat(Results_LAB_CCD, "mu_logMratio")),
              unlist(extract_stat(Results_LAB_CCD, "sigma_logMratio"))),
  NT_EXP2  = c(unlist(extract_stat(Results_LAB_NT,  "mu_logMratio")),
              unlist(extract_stat(Results_LAB_NT,  "sigma_logMratio"))),
  Variable = c('mu', 'lower mu', 'upper mu', 'sigma', 'lower sigma', 'upper sigma')
)
statGroup

# Print VR group means for quick reference
cat("VR CCD mu_Mratio:",    exp(Results_VR_CCD[[3]]$mean[Results_VR_CCD[[3]]$name == "mu_logMratio"]),   "\n")
cat("VR NT  mu_Mratio:",    exp(Results_VR_NT[[3]]$mean[Results_VR_NT[[3]]$name  == "mu_logMratio"]),    "\n")
cat("VR CCD lat mu_Mratio:",exp(Results_VR_CCD_L[[3]]$mean[Results_VR_CCD_L[[3]]$name == "mu_logMratio"]),"\n")
cat("VR NT  lat mu_Mratio:",exp(Results_VR_NT_L[[3]]$mean[Results_VR_NT_L[[3]]$name   == "mu_logMratio"]),"\n")
cat("VR CCD mon mu_Mratio:",exp(Results_VR_CCD_M[[3]]$mean[Results_VR_CCD_M[[3]]$name == "mu_logMratio"]),"\n")
cat("VR NT  mon mu_Mratio:",exp(Results_VR_NT_M[[3]]$mean[Results_VR_NT_M[[3]]$name   == "mu_logMratio"]),"\n")

# --- 6e. Posterior distribution plots (mu_logMratio) -----------
# Shared theme for HMeta-D histogram panels
hmeta_theme <- list(
  coord_cartesian(xlim = c(0, 1.5)),
  scale_y_continuous(expand = c(0, 0)),
  ggdist::theme_tidybayes(),
  theme(panel.background = element_rect(fill = 'white'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'))
)

# Experiments 1 & 2
p1 <- mcmcGroup %>%
  filter(Parameter == "mu_logMratio") %>%
  ggplot(aes(exp(value), fill = group)) +
  geom_histogram(binwidth = 0.03, alpha = 0.9, position = 'dodge') +
  scale_fill_manual(values = colpal) +
  labs(y = "Sample Count", x = expression(paste(mu["ratio"]))) +
  hmeta_theme +
  theme(legend.position = c(0.15, 0.70), legend.background = element_blank(),
        legend.title = element_blank())

# Experiment 3: binocular, lateralized, monocular panels
p1b <- make_vr_hist(c('CCD (VR)',   'NT (VR)'))
p1c <- make_vr_hist(c('CCD (VR L)', 'NT (VR L)'))
p1d <- make_vr_hist(c('CCD (VR M)', 'NT (VR M)'))

p1; p1b; p1c; p1d

# ── 7. INDIVIDUAL-LEVEL META-D' ESTIMATES -----------------------------
# Compute individual d', meta-d', and M-ratio values using make_results()
# and visualise spread across participants and presentation types.

# NOTE: make_results() is slow; gate with do == 1 to re-run

do_meta_d_d <- 0
if (do_meta_d_d == 1) {

  setwd('HMeta-d-master/R/')
  source('Function_metad_group.R')

  ccd_onl <- make_results(ONL_CCD_metaDat,   'CCD', 'ONL')
  ccd_lab <- make_results(LAB_CCD_metaDat,   'CCD', 'LAB')
  ccd_b   <- make_results(VR_CCD_metaDat,    'CCD', 'BINO')
  ccd_l   <- make_results(VR_CCD_metaDat_L,  'CCD', 'LAT')
  ccd_m   <- make_results(VR_CCD_metaDat_M,  'CCD', 'MONO')

  nt_onl  <- make_results(ONL_NT_metaDat,    'NT',  'ONL')
  nt_lab  <- make_results(LAB_NT_metaDat,    'NT',  'LAB')
  nt_b    <- make_results(VR_NT_metaDat,     'NT',  'BINO')
  nt_l    <- make_results(VR_NT_metaDat_L,   'NT',  'LAT')
  nt_m    <- make_results(VR_NT_metaDat_M,   'NT',  'MONO')

  meta_d_df <- rbind(ccd_onl, ccd_lab, ccd_b, ccd_l, ccd_m,
                      nt_onl,  nt_lab,  nt_b,  nt_l,  nt_m)

  setwd("~")
  saveRDS(meta_d_df,   'HMetaDFit/Meta_d_d/meta_d_df.rdata')

}

meta_d_df <- readRDS('HMetaDFit/Meta_d_d/meta_d_df.rdata')

# --- 7a. LM tests of group differences ---------
# Online and Lab for M-ratio, meta-d', d'
for (meas in c('Mratio', 'meta_d', 'd1')) {
  for (env in c('ONL', 'LAB', 'BINO')) {
    cat("\n", meas, env, "\n")
    print(summary(lm(mean ~ group,
                     data = meta_d_df %>%
                       filter(str_detect(name, meas),
                              !name %in% c('mu_logMratio', 'sigma_logMratio', 'mu_meta_d', 'mu_d1'),
                              type == env))))
  }
}

# --- 7b. Group-level summary plots -----
d_prime_plot_1 <- ggplot(meta_d_df %>% filter(str_detect(name, "d1"), !name %in% c('mu_d1'), type %in% c('ONL','LAB')),
                          aes(type, mean, fill = group))
meta_d_plot_1  <- ggplot(meta_d_df %>% filter(str_detect(name, "meta_d"), !name %in% c('mu_meta_d'), type %in% c('ONL','LAB')),
                          aes(type, mean, fill = group))
m_ratio_plot_1 <- ggplot(meta_d_df %>% filter(str_detect(name, "Mratio"), !name %in% c('mu_logMratio', 'sigma_logMratio'), type %in% c('ONL','LAB')),
                          aes(type, mean, fill = group))

(d_prime_plot_1 | meta_d_plot_1 | m_ratio_plot_1) &
  stat_summary(shape = 21, size = 1.2) &
  coord_cartesian(ylim = c(0, 3)) &
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 3, 0.5)) &
  scale_fill_manual(values = colpal) &
  theme_bw(base_size = 18) &
  theme(legend.position = 'none', panel.grid = element_blank(), axis.title = element_blank())

d_prime_plot <- ggplot(meta_d_df %>% filter(str_detect(name, "d1"), !str_detect(name, "mu_d1"), type %in% c('BINO','MONO','LAT')),
                        aes(type, abs(mean), fill = group))
meta_d_plot  <- ggplot(meta_d_df %>% filter(str_detect(name, "meta_d"), !str_detect(name, "mu_meta_d"), type %in% c('BINO','MONO','LAT')),
                        aes(type, abs(mean), fill = group))
m_ratio_plot <- ggplot(meta_d_df %>% filter(str_detect(name, "Mratio"), !str_detect(name, "mu_logMratio"), !str_detect(name, "sigma_logMratio"), type %in% c('BINO','MONO','LAT')),
                        aes(type, abs(mean), fill = group))

(d_prime_plot | meta_d_plot | m_ratio_plot) &
  stat_summary(geom = 'line', aes(group = group)) &
  stat_summary(shape = 21, size = 1.2) &
  coord_cartesian(ylim = c(0, 2)) &
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 2, 0.5)) &
  scale_fill_manual(values = colpal_vr[1:2]) &
  theme_bw(base_size = 18) &
  theme(legend.position = 'none', panel.grid = element_blank(), axis.title = element_blank())

# ── 8. MODEL RECOVERY & POSTERIOR PREDICTIVE CHECK ---------------------
# Validate the HMeta-D model by recovering the fitted response counts from
# the posterior and correlating them with the empirical counts.

## NOTE: outcomes are stochastic. Expect some deviation from published
## results if running again.

# --- 8a. Run recovery -----

do <- 0
if(do==1){

  setwd('HMeta-d-master/R/')
  source('Function_metad_group.R')

  fit_CCD_lab_rec <- run_recovery(LAB_CCD_metaDat, 'CCD', 'LAB')
  fit_NT_lab_rec  <- run_recovery(LAB_NT_metaDat,  'NT',  'LAB')
  fit_CCD_onl_rec <- run_recovery(ONL_CCD_metaDat, 'CCD', 'ONL')
  fit_NT_onl_rec  <- run_recovery(ONL_NT_metaDat,  'NT',  'ONL')
  fit_CCD_bin_rec <- run_recovery(VR_CCD_metaDat,  'CCD', 'BINO')
  fit_NT_bin_rec  <- run_recovery(VR_NT_metaDat,   'NT',  'BINO')
  fit_CCD_lat_rec <- run_recovery(VR_CCD_metaDat_L,'CCD', 'LAT')
  fit_NT_lat_rec  <- run_recovery(VR_NT_metaDat_L, 'NT',  'LAT')
  fit_CCD_mon_rec <- run_recovery(VR_CCD_metaDat_M,'CCD', 'MONO')
  fit_NT_mon_rec  <- run_recovery(VR_NT_metaDat_M, 'NT',  'MONO')

  recovered_dfs <- rbind(fit_CCD_lab_rec[[1]], fit_NT_lab_rec[[1]],
                        fit_CCD_onl_rec[[1]], fit_NT_onl_rec[[1]],
                        fit_CCD_bin_rec[[1]], fit_NT_bin_rec[[1]]) %>%
  mutate(type = case_when(type == 'ONL'  ~ 'EXP1',
                          type == 'LAB'  ~ 'EXP2',
                          type == 'BINO' ~ 'EXP3',
                          .default = type))

  setwd("~")
  saveRDS(recovered_dfs,   'HMetaDFit/Recovery/recovered_dfs.rdata')

}

recovered_dfs <- readRDS('HMetaDFit/Recovery/recovered_dfs.rdata')

# Correlation of real vs recovered M-ratio per participant
ggplot(recovered_dfs %>% filter(str_detect(name,'Mratio'), !str_detect(name,'log')),
       aes(mean_real, mean_rec)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', se = FALSE) +
  ggpubr::stat_cor(method = 'spearman', label.x.npc = 0.25, label.y.npc = 0,
                   colour = 'firebrick', size = 5) +
  facet_wrap(group ~ type, scales = 'free', nrow = 2) +
  labs(x = 'Real', y = 'Recovered') +
  theme_minimal(base_size = 18) + theme(panel.grid = element_blank())

# --- 8b. Posterior predictive check -----------
# Compare recovered and real response count distributions across confidence bins
post_pred_check <- data.frame(
  conf  = rep(c(rep(c(6:1,-1:-6), 4), rep(c(2:1,-1:-2), 2)), 2),
  cond  = rep(c(rep('EXP2',24), rep('EXP1',24), rep('EXP3',8)), 2),
  type  = c(rep('Rec',56), rep('Real',56)),
  group = rep(c(rep(c(rep('CCD',12), rep('NT',12)), 2), c(rep('CCD',4), rep('NT',4))), 2),
  s1    = c(rowSums(fit_CCD_lab_rec[[2]][[1]]), rowSums(fit_NT_lab_rec[[2]][[1]]),
            rowSums(fit_CCD_onl_rec[[2]][[1]]), rowSums(fit_NT_onl_rec[[2]][[1]]),
            rowSums(fit_CCD_bin_rec[[2]][[1]]), rowSums(fit_NT_bin_rec[[2]][[1]]),
            rowSums(LAB_CCD_metaDat[[1]]),      rowSums(LAB_NT_metaDat[[1]]),
            rowSums(ONL_CCD_metaDat[[1]]),      rowSums(ONL_NT_metaDat[[1]]),
            rowSums(VR_CCD_metaDat[[1]]),       rowSums(VR_NT_metaDat[[1]])),
  s2    = c(rowSums(fit_CCD_lab_rec[[2]][[2]]), rowSums(fit_NT_lab_rec[[2]][[2]]),
            rowSums(fit_CCD_onl_rec[[2]][[2]]), rowSums(fit_NT_onl_rec[[2]][[2]]),
            rowSums(fit_CCD_bin_rec[[2]][[2]]), rowSums(fit_NT_bin_rec[[2]][[2]]),
            rowSums(LAB_CCD_metaDat[[2]]),      rowSums(LAB_NT_metaDat[[2]]),
            rowSums(ONL_CCD_metaDat[[2]]),      rowSums(ONL_NT_metaDat[[2]]),
            rowSums(VR_CCD_metaDat[[2]]),       rowSums(VR_NT_metaDat[[2]]))
) %>%
  pivot_longer(s1:s2, names_to = 'Dec', values_to = 'Count') %>%
  group_by(type, Dec, cond, group) %>%
  mutate(Count = Count / sum(Count))

# Compute per-bin Pearson r between real and recovered counts for each experiment
bin_correlations <- bind_rows(
  bind_rows(fit_to_count(fit_CCD_onl_rec[[2]][[1]],"Rec","EXP1","s1"),
            fit_to_count(ONL_CCD_metaDat[[1]],      "Real","EXP1","s1"),
            fit_to_count(fit_NT_onl_rec[[2]][[2]],  "Rec","EXP1","s2"),
            fit_to_count(ONL_NT_metaDat[[2]],        "Real","EXP1","s2")) %>% mutate(exp='EXP1'),
  bind_rows(fit_to_count(fit_CCD_lab_rec[[2]][[1]],"Rec","EXP2","s1"),
            fit_to_count(LAB_CCD_metaDat[[1]],      "Real","EXP2","s1"),
            fit_to_count(fit_NT_lab_rec[[2]][[2]],  "Rec","EXP2","s2"),
            fit_to_count(LAB_NT_metaDat[[2]],        "Real","EXP2","s2")) %>% mutate(exp='EXP2'),
  bind_rows(fit_to_count(fit_CCD_bin_rec[[2]][[1]],"Rec","EXP3","s1"),
            fit_to_count(VR_CCD_metaDat[[1]],       "Real","EXP3","s1"),
            fit_to_count(fit_NT_bin_rec[[2]][[2]],  "Rec","EXP3","s2"),
            fit_to_count(VR_NT_metaDat[[2]],         "Real","EXP3","s2")) %>% mutate(exp='EXP3')
) %>%
  pivot_wider(names_from = Type, values_from = Count) %>%
  group_by(Bin, exp) %>%
  summarise(test_results = list(broom::tidy(cor.test(Real, Rec, method="pearson"))),
            .groups = 'drop') %>%
  unnest(test_results) %>%
  dplyr::select(Bin, exp, pearson_r = estimate, p_value = p.value,
                conf_low = conf.low, conf_high = conf.high) %>%
  mutate(conf_level = case_when(
    exp == 'EXP3' & Bin == 1 ~  2, exp == 'EXP3' & Bin == 2 ~  1,
    exp == 'EXP3' & Bin == 3 ~ -1, exp == 'EXP3' & Bin == 4 ~ -2,
    exp %in% c('EXP1','EXP2') & Bin %in% 1:6  ~  7 - Bin,
    exp %in% c('EXP1','EXP2') & Bin %in% 7:12 ~  6 - Bin,
    .default = as.numeric(Bin)),
    conf_level = as.numeric(conf_level))

bin_correlations %>%
  filter(exp=='EXP3')

post_pred_check %>%
  pivot_wider(names_from = type, values_from = Count) %>%

  ggplot(aes(x = Real, y = Rec, color = factor(conf))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + # Identity line
  facet_wrap(~cond) +
  labs(x = "Real",
       y = "Rec.")+
  theme_bw(base_size = 24)+
  theme(legend.position = 'top',
          legend.background = element_blank(),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          strip.background.x = element_blank())

# Plot raw distribution comparison + per-bin r (stacked)
post_pred_check_plot <- ggplot(post_pred_check, aes(conf, Count, fill = type)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('#2D3142','#F6AE2D')) +
  scale_x_continuous(breaks = c(6,1,-1,-6)) +
  scale_y_continuous(breaks = seq(0.1,0.5,0.1), expand = c(0,0)) +
  facet_wrap(group ~ cond, scales = 'free') +
  labs(y = 'Normalised Count') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none', panel.grid = element_blank(),
        strip.background.x = element_blank(), axis.title.x = element_blank())

cor_plot_post_pred <- ggplot(bin_correlations, aes(conf_level, pearson_r)) +
  geom_col(fill = "black") +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Confidence Bin (1=low, 6=high)", y = "Pearson r") +
  facet_wrap(~exp, scales = 'free') +
  scale_x_continuous(breaks = c(6,1,-1,-6)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), strip.background.x = element_blank(),
        strip.text.x = element_blank())

(post_pred_check_plot / cor_plot_post_pred) & plot_layout(heights = c(5,1))

# ── 9. PERMUTATION TESTS ON META-D' GROUP DIFFERENCES -------------------
# We use a permutation procedure to test whether the observed NT − CCD difference
# in mu_logMratio (and d', meta-d') exceeds the null distribution
# generated by label-shuffling.

# --- 9a. Extract observed group differences ----
# Helper: pull exp(mu_logMratio) from a make_results() output dataframe
pull_mu <- function(df, param = 'mu_logMratio') df %>% filter(name == param) %>% mutate(mean = exp(mean))
pull_meta_d <- function(df, param = 'mu_meta_d') df %>% filter(name == param) %>% mutate(mean = abs(mean))
pull_d1 <- function(df, param = 'mu_d1') df %>% filter(name == param) %>% mutate(mean = abs(mean))

mu_log_exp1_NT  <- pull_mu(Results_ONL_NT$Fit);      mu_log_exp1_CCD  <- pull_mu(Results_ONL_CCD$Fit)
mu_log_exp2_NT  <- pull_mu(Results_LAB_NT$Fit);      mu_log_exp2_CCD  <- pull_mu(Results_LAB_CCD$Fit)
mu_log_exp3_NT   <- pull_mu(Results_VR_NT$Fit);       mu_log_exp3_CCD   <- pull_mu(Results_VR_CCD$Fit)
mu_log_exp3_NT_L <- pull_mu(Results_VR_NT_L$Fit);     mu_log_exp3_CCD_L <- pull_mu(Results_VR_CCD_L$Fit)
mu_log_exp3_NT_M <- pull_mu(Results_VR_NT_M$Fit);     mu_log_exp3_CCD_M <- pull_mu(Results_VR_CCD_M$Fit)

meta_d_log_exp1_NT   <- pull_meta_d(nt_onl);meta_d_log_exp1_CCD   <- pull_meta_d(ccd_onl);
meta_d_log_exp2_NT   <- pull_meta_d(nt_lab);meta_d_log_exp2_CCD   <- pull_meta_d(ccd_lab);
meta_d_log_exp3_NT   <- pull_meta_d(nt_b);  meta_d_log_exp3_CCD   <- pull_meta_d(ccd_b);
meta_d_log_exp3_NT_L <- pull_meta_d(nt_l);  meta_d_log_exp3_CCD_L <- pull_meta_d(ccd_l);
meta_d_log_exp3_NT_M <- pull_meta_d(nt_m);  meta_d_log_exp3_CCD_M <- pull_meta_d(ccd_m);

d1_log_exp1_NT   <- pull_d1(nt_onl);d1_log_exp1_CCD   <- pull_d1(ccd_onl);
d1_log_exp2_NT   <- pull_d1(nt_lab);d1_log_exp2_CCD   <- pull_d1(ccd_lab);
d1_log_exp3_NT   <- pull_d1(nt_b);  d1_log_exp3_CCD   <- pull_d1(ccd_b);
d1_log_exp3_NT_L <- pull_d1(nt_l);  d1_log_exp3_CCD_L <- pull_d1(ccd_l);
d1_log_exp3_NT_M <- pull_d1(nt_m);  d1_log_exp3_CCD_M <- pull_d1(ccd_m);

# Compute NT − CCD differences for M-ratio, d', and meta-d' in each context
diffExp1         <- as.numeric(mu_log_exp1_NT[2]   - mu_log_exp1_CCD[2])
diffExp2         <- as.numeric(mu_log_exp2_NT[2]   - mu_log_exp2_CCD[2])
diffExp3         <- as.numeric(mu_log_exp3_NT[2]   - mu_log_exp3_CCD[2])
diffExp3_L       <- as.numeric(mu_log_exp3_NT_L[2] - mu_log_exp3_CCD_L[2])
diffExp3_M       <- as.numeric(mu_log_exp3_NT_M[2] - mu_log_exp3_CCD_M[2])

diffExp1_d1         <- as.numeric(d1_log_exp1_NT[2]   - d1_log_exp1_CCD[2])
diffExp2_d1         <- as.numeric(d1_log_exp2_NT[2]   - d1_log_exp2_CCD[2])
diffExp3_d1         <- as.numeric(d1_log_exp3_NT[2]   - d1_log_exp3_CCD[2])
diffExp3_L_d1       <- as.numeric(d1_log_exp3_NT_L[2] - d1_log_exp3_CCD_L[2])
diffExp3_M_d1       <- as.numeric(d1_log_exp3_NT_M[2] - d1_log_exp3_CCD_M[2])

diffExp1_metad         <- as.numeric(meta_d_log_exp1_NT[2]   - meta_d_log_exp1_CCD[2])
diffExp2_metad         <- as.numeric(meta_d_log_exp2_NT[2]   - meta_d_log_exp2_CCD[2])
diffExp3_metad         <- as.numeric(meta_d_log_exp3_NT[2]   - meta_d_log_exp3_CCD[2])
diffExp3_L_metad       <- as.numeric(meta_d_log_exp3_NT_L[2] - meta_d_log_exp3_CCD_L[2])
diffExp3_M_metad       <- as.numeric(meta_d_log_exp3_NT_M[2] - meta_d_log_exp3_CCD_M[2])

cat("mu difference Lab:",  diffExp1,  "\n")
cat("mu difference Onl:",  diffExp2,  "\n")
cat("mu difference VR:",   diffExp3,   "\n")
cat("mu difference VR-L:", diffExp3_L, "\n")
cat("mu difference VR-M:", diffExp3_M, "\n")

# --- 9b. Prepare data for permutation ----
data_for_permute_exp1_CCD      <- checkBothRDK %>% filter(group=='CCD (Online)') %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)
data_for_permute_exp1_NT       <- checkBothRDK %>% filter(group=='NT (Online)')  %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)

data_for_permute_exp2_CCD      <- checkBothRDK %>% filter(group=='CCD (MRI)')    %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)
data_for_permute_exp2_NT       <- checkBothRDK %>% filter(group=='NT (MRI)')     %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)

data_for_permute_exp3_CCD      <- data_vr %>% filter(BinocularTrial   == 1, Group=='CCD')
data_for_permute_exp3_NT       <- data_vr %>% filter(BinocularTrial   == 1, Group=='NT')
data_for_permute_exp3_CCD_L    <- data_vr %>% filter(LateralizedTrial == 1, Group=='CCD')
data_for_permute_exp3_NT_L     <- data_vr %>% filter(LateralizedTrial == 1, Group=='NT')
data_for_permute_exp3_CCD_M    <- data_vr %>% filter(MonocularTrial   == 1, Group=='CCD')
data_for_permute_exp3_NT_M     <- data_vr %>% filter(MonocularTrial   == 1, Group=='NT')

# --- 9c. Run permutation (slow; only re-run when needed) ----

## NOTE: outcomes are stochastic. Expect some deviation from published
## results if running again.
do_perm <- 0
if(do_perm==1)

  setwd('HMeta-d-master/R/')
  source('Function_metad_group.R')
  nreps <- 500

  onl_CCDvsNTp      <- permute_metad_group(data_for_permute_exp1_CCD, data_for_permute_exp1_NT, cores = 5, nreps = nreps);
  onl_CCDvsNTp_mrat <- permute_metad_group_mrat(data_for_permute_exp1_CCD, data_for_permute_exp1_NT, cores = 5, nreps = nreps);

  lab_CCDvsNTp      <- permute_metad_group(data_for_permute_exp2_CCD, data_for_permute_exp2_NT, cores = 10, nreps = nreps);
  lab_CCDvsNTp_mrat <- permute_metad_group_mrat(data_for_permute_exp2_CCD, data_for_permute_exp2_NT, cores = 10, nreps = nreps);

  vr_CCDvsNTp  <- permute_metad_group_vr(data_for_permute_exp3_CCD,   data_for_permute_exp3_NT,cores = 10, nreps = nreps);
  vr_CCDvsNTp_L<- permute_metad_group_vr(data_for_permute_exp3_CCD_L, data_for_permute_exp3_NT_L,cores = 10, nreps = nreps);
  vr_CCDvsNTp_M<- permute_metad_group_vr(data_for_permute_exp3_CCD_M, data_for_permute_exp3_NT_M,cores = 10, nreps = nreps);

  saveRDS(onl_CCDvsNTp,  'HMetaDFit/CCDvsNT_ONL.rdata')
  saveRDS(onl_CCDvsNTp,  'HMetaDFit/CCDvsNT_ONL_mrat.rdata')
  saveRDS(lab_CCDvsNTp,  'HMetaDFit/CCDvsNT_LAB.rdata')
  saveRDS(lab_CCDvsNTp_mrat,  'HMetaDFit/CCDvsNT_LAB_mrat.rdata')
  saveRDS(vr_CCDvsNTp,  'HMetaDFit/CCDvsNT_VR.rdata')
  saveRDS(vr_CCDvsNTp_L,'HMetaDFit/CCDvsNT_VR_L.rdata')
  saveRDS(vr_CCDvsNTp_M,'HMetaDFit/CCDvsNT_VR_M.rdata')

  setwd('~')

}

# --- 9d. Load permutation null distributions ----
onl_CCDvsNTp      <- readRDS('HMetaDFit/Permutation/CCDvsNT_ONL.rdata')
onl_CCDvsNTp_mrat <- readRDS('HMetaDFit/Permutation/CCDvsNT_ONL_mrat.rdata')
lab_CCDvsNTp      <- readRDS('HMetaDFit/Permutation/CCDvsNT_LAB.rdata')
lab_CCDvsNTp_mrat <- readRDS('HMetaDFit/Permutation/CCDvsNT_LAB_mrat.rdata')
vr_CCDvsNTp       <- readRDS('HMetaDFit/Permutation/CCDvsNT_VR.rdata')
vr_CCDvsNTp_L     <- readRDS('HMetaDFit/Permutation/CCDvsNT_VR_L.rdata')
vr_CCDvsNTp_M     <- readRDS('HMetaDFit/Permutation/CCDvsNT_VR_M.rdata')

# --- 9e. Effect sizes and permutation p-values ------

# perm_effect_stats() returns the p-value and U3 superiority index
cat("\n--- Online (Exp 1) ---\n")
perm_effect_stats(round(diffExp1,         2), onl_CCDvsNTp_mrat[,4])
perm_effect_stats(round(diffExp1_d1,      2), onl_CCDvsNTp[,7])
perm_effect_stats(round(diffExp1_metad,   2), onl_CCDvsNTp[,10])

cat("\n--- Lab (Exp 2) ---\n")
perm_effect_stats(round(diffExp2,         2), lab_CCDvsNTp_mrat[,4])
perm_effect_stats(round(diffExp2_d1,      2), lab_CCDvsNTp[,7])
perm_effect_stats(round(diffExp2_metad,   2), lab_CCDvsNTp[,10])

cat("\n--- VR binocular (Exp 3) ---\n")
perm_effect_stats(round(diffExp3,         2), vr_CCDvsNTp_mrat[,4])
perm_effect_stats(round(diffExp3_d1,      2), vr_CCDvsNTp[,7])
perm_effect_stats(round(diffExp3_meta_d,  2), vr_CCDvsNTp[,10])

cat("\n--- VR lateralized (Exp 3) ---\n")
perm_effect_stats(round(diffExp3_L,       2), vr_CCDvsNTp_mrat_L[,4])
perm_effect_stats(round(diffExp3_L_d1,    2), vr_CCDvsNTp_L[,7])
perm_effect_stats(round(diffExp3_L_metad, 2), vr_CCDvsNTp_L[,10])

cat("\n--- VR monocular (Exp 3) ---\n")
perm_effect_stats(round(diffExp3_M,       2), vr_CCDvsNTp_mrat_M[,4])
perm_effect_stats(round(diffExp3_M_d1,    2), vr_CCDvsNTp_M[,7])
perm_effect_stats(round(diffExp3_M_metad, 2), vr_CCDvsNTp_M[,10])

# --- 9f. Permutation histogram plots -------
# One panel per experiment, x-axis = null M-ratio difference, vertical line = observed
p_p1 <- onl_CCDvsNTp %>% as.data.frame() %>% rename(Diff = 4) %>%
  ggplot(aes(Diff)) + geom_histogram(binwidth = 0.01) +
  labs(x = expression(paste("Randomised ",mu['ratio'])), y = 'Sample Count') +
  geom_label(x = diffExp1 + 0.07, y = 10, label = round(diffExp1,2), fill='#81D6E3', colour='white') +
  geom_vline(xintercept = diffExp1, color = '#81D6E3', size = 1.5) +
  theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_text(size=24), axis.title = element_text(size=24))

p_p2 <- lab_CCDvsNTp %>% as.data.frame() %>% rename(Diff = 4) %>%
  ggplot(aes(Diff)) + geom_histogram(binwidth = 0.01) +
  labs(x = expression(paste("Randomised ",mu['ratio'])), y = 'Sample Count') +
  geom_label(x = diffExp2 + 0.07, y = 10, label = round(diffExp2,2), fill='#475B5A', colour='white') +
  geom_vline(xintercept = diffExp2, color = '#475B5A', size = 1.5) +
  theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_text(size=24), axis.title = element_text(size=24))

p_p_b <- make_vr_perm_plot(vr_CCDvsNTp,   diffExp3)
p_p_l <- make_vr_perm_plot(vr_CCDvsNTp_L, diffExp3_L)
p_p_m <- make_vr_perm_plot(vr_CCDvsNTp_M, diffExp3_M)

(p_p1 | p_p2)
(p_p_b | p_p_l | p_p_m)

# ── 10. SUPPLEMENT: REACTION TIME ANALYSIS -----------------------------
# Supplementary: RT distributions and LME tests by coherence and correctness.
# Decision trials (conf == 0) and confidence trials (conf != 0) treated separately.

library(ggridges)
Experiment_1_both <- rbind(Experiment_1_CCD, Experiment_1_NT)
Experiment_2_both <- rbind(Experiment_2_CCD, Experiment_2_NT)

## Non Confidence Trials -----

### Split by coherence ----

Experiment_1_both %>%
  filter(type == 'main', kmed != 'kmed', ID != '5f20a850b6d5ed3d186bc3c3', conf == 0) %>%
  group_by(group, kmed, ID) %>%
  summarise(Subj_Acc = mean(correct)*100, Subj_RT = mean(RT, na.rm=TRUE), .groups='drop') %>%
  group_by(group, kmed) %>%
  summarise(Mean_Acc=mean(Subj_Acc), SD_Acc=sd(Subj_Acc),
            Mean_RT=mean(Subj_RT),   SD_RT=sd(Subj_RT),
            SE_RT=SD_RT/sqrt(n()),   n=n(), .groups='drop')

Experiment_2_both %>%
  filter(type == 'main', kmed != 'kmed', conf == 0) %>%
  group_by(group, kmed, ID) %>%
  summarise(Subj_Acc = mean(correct)*100, Subj_RT = mean(RT, na.rm=TRUE), .groups='drop') %>%
  group_by(group, kmed) %>%
  summarise(Mean_Acc=mean(Subj_Acc), SD_Acc=sd(Subj_Acc),
            Mean_RT=mean(Subj_RT),   SD_RT=sd(Subj_RT),
            SE_RT=SD_RT/sqrt(n()),   n=n(), .groups='drop')

run_rt_lmer(Experiment_1_both, 'Computer', conf_filter = 0)  # Exp 1: decision RT × coherence
run_rt_lmer(Experiment_2_both, 'MRI',      conf_filter = 0)  # Exp 2: decision RT × coherence

### Split by correctness ----

Experiment_1_both %>%
  filter(type == 'main', kmed != 'kmed', ID != '5f20a850b6d5ed3d186bc3c3', conf == 0) %>%
  group_by(group, correct, ID) %>%
  summarise(Subj_Acc = mean(correct)*100, Subj_RT = mean(RT, na.rm=TRUE), .groups='drop') %>%
  group_by(group, correct) %>%
  summarise(Mean_Acc=mean(Subj_Acc), SD_Acc=sd(Subj_Acc),
            Mean_RT=mean(Subj_RT),   SD_RT=sd(Subj_RT),
            SE_RT=SD_RT/sqrt(n()),   n=n(), .groups='drop')

Experiment_2_both %>%
  filter(type == 'main', kmed != 'kmed', conf == 0) %>%
  group_by(group, correct, ID) %>%
  summarise(Subj_Acc = mean(correct)*100, Subj_RT = mean(RT, na.rm=TRUE), .groups='drop') %>%
  group_by(group, correct) %>%
  summarise(Mean_Acc=mean(Subj_Acc), SD_Acc=sd(Subj_Acc),
            Mean_RT=mean(Subj_RT),   SD_RT=sd(Subj_RT),
            SE_RT=SD_RT/sqrt(n()),   n=n(), .groups='drop')

lmerTest::lmer(log(RT) ~ CCD * correct + (1|ID),
               data = Experiment_1_both %>% filter(type=='main', conf==0, kmed!='kmed', environment=='Computer')) %>% summary()
lmerTest::lmer(log(RT) ~ CCD * correct + (1|ID),
               data = Experiment_2_both %>% filter(type=='main', conf==0, kmed!='kmed', environment=='MRI'))      %>% summary()

# VR RT models by presentation type
data_vr %>%
  filter(TrialType == 'Main') %>%
  group_by(Group, PresentationType, HighCoherenceTrialLabel, PID) %>%
  summarise(Subj_Acc = mean(Correct)*100, Subj_RT = mean(ReactionTime, na.rm=TRUE), .groups='drop') %>%
  group_by(Group, PresentationType, HighCoherenceTrialLabel) %>%
  summarise(Mean_Acc=mean(Subj_Acc), SD_Acc=sd(Subj_Acc),
            Mean_RT=mean(Subj_RT),   SD_RT=sd(Subj_RT),
            SE_RT=SD_RT/sqrt(n()),   n=n(), .groups='drop')

for (ptype in c('Binocular','Lateralized','Monocular')) {
  cat("\n=== RT VR", ptype, "===\n")
  lmerTest::lmer(log(ReactionTime) ~ CCD * HighCoherenceTrialLabel + (1|PID),
                 data = data_vr %>% filter(TrialType=='Main', PresentationType==ptype,
                                           HighCoherenceTrialLabel %in% c('x0.5','x2'))) %>% summary() %>% print()
  lmerTest::lmer(log(ReactionTime) ~ CCD * Correct + (1|PID),
                 data = data_vr %>% filter(TrialType=='Main', PresentationType==ptype,
                                           HighCoherenceTrialLabel %in% c('x0.5','x2'))) %>% summary() %>% print()
}

# ── 11. SUPPLEMENT: CALIBRATION ASSESSMENT ---------------------------------
# Show that k_med staircase converges similarly in both groups and environments.

## Non-VR Based ----

ID_baseK <- checkBothRDK %>%
  filter(type == 'calibration', Trial > 100) %>%
  group_by(ID, CCD, environment) %>%
  summarise(baseK = median(coherence))

allNonVRK <- ID_baseK %>%
  group_by(CCD, environment) %>%
  summarise(Coherence = mean(baseK),
            CoherenceSD = sd(baseK))
allNonVRK
# Test whether base k varied across computer or MRI versions of the same task
t.test(x = (ID_baseK %>% filter(CCD == 1, environment == 'Computer'))$baseK,
       y = (ID_baseK %>% filter(CCD == 0, environment == 'Computer'))$baseK)
t.test(x = (ID_baseK %>% filter(CCD == 1, environment == 'MRI'))$baseK,
       y = (ID_baseK %>% filter(CCD == 0, environment == 'MRI'))$baseK)

# Base K between conditions ------
coh_tog <-   rbind(
  ID_baseK,
  data_vr %>%
        filter(TrialType == 'Main') %>%
        group_by(PID, Group) %>%
        summarise(Coherence = mean(BaseCoherence)) %>%
        mutate(environment = 'VR',
               Group = ifelse(Group=='CCD', 1, 0),
               PID = as.character(PID)) %>%
        dplyr::select(ID = PID, CCD = Group, environment, baseK = Coherence)
  )

aov_tog_ccd <- aov(baseK ~ environment, coh_tog %>% filter(CCD==1))
summary(aov_tog_ccd)
TukeyHSD(aov_tog_ccd)

aov_tog_nt <- aov(baseK ~ environment, coh_tog %>% filter(CCD==0))
summary(aov_tog_nt)
TukeyHSD(aov_tog_nt)

## VR Based -----
allVRK <- data_vr %>%
  filter(TrialType == 'Main') %>%
  group_by(CCD, PresentationType) %>%
  summarise(Coherence = mean(BaseCoherence),
            CoherenceSD = sd(BaseCoherence))

binoK <- data_vr %>%
  filter(BinocularTrial == 1, TrialType == 'Main') %>%
  group_by(PID, CCD) %>%
  summarise(Coherence = mean(BaseCoherence))

monoK <- data_vr %>%
  filter(MonocularTrial == 1, TrialType == 'Main') %>%
  group_by(PID, CCD) %>%
  summarise(Coherence = mean(BaseCoherence))

latK <- data_vr %>%
  filter(LateralizedTrial == 1, TrialType == 'Main') %>%
  group_by(PID, CCD) %>%
  summarise(Coherence = mean(BaseCoherence))

# Test whether base k varied across presentations
t.test((binoK %>% filter(CCD == 1))$Coherence, y = (binoK %>% filter(CCD == 0))$Coherence)
t.test((monoK %>% filter(CCD == 1))$Coherence, y = (monoK %>% filter(CCD == 0))$Coherence)
t.test((latK %>% filter(CCD == 1))$Coherence, y = (latK %>% filter(CCD == 0))$Coherence)

#ANOVA
coh_aov_ccd <- lmerTest::lmer(Coherence ~ PresentationType + (1|PID),
    data =
      data_vr %>%
      filter(TrialType == 'Main', Group == 'CCD') %>%
      group_by(PID, PresentationType) %>%
      mutate(PresentationType = factor(PresentationType, levels = c( 'Monocular', 'Binocular', 'Lateralized'))) %>%
      summarise(Coherence = mean(BaseCoherence)))
anova(coh_aov_ccd)
summary(coh_aov_ccd)

coh_aov_nt <- lmerTest::lmer(Coherence ~ PresentationType + (1|PID),
    data =
      data_vr %>%
      filter(TrialType == 'Main', Group == 'NT') %>%
      group_by(PID, PresentationType) %>%
      summarise(Coherence = mean(BaseCoherence)))
anova(coh_aov_nt)
summary(coh_aov_nt)

coh_aov_comb <- lmerTest::lmer(Coherence ~ PresentationType * Group + (1|PID),
    data =
      data_vr %>%
      filter(TrialType == 'Main') %>%
      group_by(PID, PresentationType, Group) %>%
      summarise(Coherence = mean(BaseCoherence)))
anova(coh_aov_comb)
summary(coh_aov_comb)


# Combine all experiments into one long format
check_coh_RDK <- rbind(
  Experiment_1_CCD, Experiment_2_CCD, Experiment_1_NT, Experiment_2_NT
) %>%
  mutate(prestype = NA, vr = 0) %>%
  dplyr::select(ID, Trial, coherence, group, prestype, type, vr) %>%
  rbind(
    data_vr %>%
      rename(ID = PID, group = Group, coherence = ActiveCoherence,
             prestype = PresentationType, type = TrialType) %>%
      mutate(vr = 1) %>%
      dplyr::select(ID, Trial, coherence, group, prestype, type, vr)
  )

# Coherence trajectory plot: Experiments 1 & 2
CCDcheck2 <-
  ggplot(check_coh_RDK %>% filter(type == 'calibration', vr == 0)) +
    stat_summary(aes(Trial, coherence, fill = group), geom='ribbon', alpha=0.1, colour=NA) +
    stat_summary(aes(Trial, coherence, colour = group), geom='line') +
    scale_color_manual(values = colpal) + scale_fill_manual(values = colpal) +
    labs(x = 'Trial', y = 'Coherence') +
    coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 118)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_minimal(base_size = 20) +
    theme(legend.position = 'top', panel.grid.minor = element_blank(), legend.title = element_blank())

bino_check <- make_vr_cal_plot('Binocular',  'Bino.')
lat_check  <- make_vr_cal_plot('Lateralized','Lat.')
mono_check <- make_vr_cal_plot('Monocular',  'Mono.')

coh_vr_check_full <- ggarrange(bino_check, lat_check, mono_check, nrow = 1)

CCDcheck2; coh_vr_check_full

# Correlation: does mean confidence predict k_med? (supplementary validity check)
check_coh_lm <- checkBothRDK %>%
  group_by(ID) %>%
  filter(type == 'main', conf != 0) %>%
  mutate(kmed_av = ifelse(kmed == 'kmed x 0.5', coherence*2, coherence/2),
         conf    = mean(conf)) %>%
  dplyr::select(kmed_av, conf, ID, group) %>%
  distinct()

ggplot(check_coh_lm, aes(conf, kmed_av, colour = group)) +
  geom_point() + geom_smooth(method = 'lm') +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_color_manual(values = colpal) +
  stat_cor(show.legend = FALSE, method = 'spearman') +
  stat_cor(show.legend = FALSE, aes(colour = NULL), data = check_coh_lm,
           label.y.npc = 0.50, method = 'spearman') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'left', legend.title = element_blank())


# ── 12. SUPPLEMENT: GROUP BEHAVIOUR FIGURES --------------------------

# Experiments 1 & 2: beeswarm accuracy and confidence plots with significance bars
tidy_acc <- tidyplot(
  av_cor %>% mutate(
    kmed  = ifelse(kmed == 'x0.5', 'x0.5', 'x2.0'),
    group = dplyr::recode(group, 'CCD (MRI)'='CCD\n(MRI)', 'NT (Online)'='NT\n(Online)',
                                  'CCD (Online)'='CCD\n(Online)', 'NT (MRI)'='NT\n(MRI)')),
  x = group, y = av_cor, colour = kmed) %>%
  add_data_points_beeswarm(alpha=0.1, cex=1) %>% add_sd_errorbar() %>% add_mean_dot(size=2) %>%
  adjust_size(width=75, height=50) %>% adjust_x_axis_title(title='') %>%
  adjust_y_axis_title(title='p(Correct)') %>% adjust_font(fontsize=15) %>%
  adjust_legend_title(title='Coh.') %>% adjust_legend_position(position='top') +
  stat_compare_means(paired=TRUE, hide.ns=TRUE, size=5, show.legend=FALSE,
                     label='p.signif', label.y=1.04) +
  coord_cartesian(clip="off")

tidy_conf <- tidyplot(
  av_co %>% mutate(
    kmed  = ifelse(kmed == 'x0.5', 'x0.5', 'x2.0'),
    group = dplyr::recode(group, 'CCD (MRI)'='CCD\n(MRI)', 'NT (Online)'='NT\n(Online)',
                                  'CCD (Online)'='CCD\n(Online)', 'NT (MRI)'='NT\n(MRI)'),
    ID    = factor(ID)),
  x = group, y = av_co, colour = kmed) %>%
  add_data_points_beeswarm(alpha=0.1, cex=1) %>% add_sd_errorbar() %>% add_mean_dot(size=2) %>%
  adjust_size(width=75, height=50) %>% adjust_x_axis_title(title='') %>%
  adjust_y_axis_title(title='Confidence') %>% adjust_font(fontsize=15) %>%
  adjust_legend_title(title='Coh.') %>% adjust_legend_position(position='top') +
  stat_compare_means(paired=TRUE, hide.ns=TRUE, size=5, show.legend=FALSE,
                     label='p.signif', label.y=102) +
  coord_cartesian(clip="off")

tidy_acc | tidy_conf

tidyplot_exp3_cor <- av_cor_vr %>%
  mutate(kmed = as.character(kmed),
         Pres = dplyr::recode(PresentationType, 'Monocular'='Mono.', 'Lateralized'='Lat.', 'Binocular'='Bino.'),
         ID   = factor(ID))

tidyplot_exp3_co <- av_co_vr %>%
  mutate(kmed = as.character(kmed),
         Pres = dplyr::recode(PresentationType, 'Monocular'='Mono.', 'Lateralized'='Lat.', 'Binocular'='Bino.'),
         ID   = factor(ID))

Cor  <- make_exp3_plot(tidyplot_exp3_cor, 'av_cor', 'p(Correct)',  1.04)
Conf <- make_exp3_plot(tidyplot_exp3_co,  'av_co',  'Confidence',  0.75)

Cor / Conf

# ── 13. SUPPLEMENT: INTELLIGENCE CORRELATION -------------------------------
# ICAR scores from online platforms linked to individual meta-d' estimates to
# check whether metacognitive efficiency correlates with fluid intelligence.

ID_no_conf <- checkBothRDK %>%
  group_by(ID) %>%
  filter(all(conf == 0)) %>%
  pull(ID) %>%
  unique()

ID_nt_mrat <- checkBothRDK %>%
  filter(group == 'NT (Online)',
         ID != ID_no_conf) %>%
  dplyr::select(ID) %>%
  unique()

ICAR_NT <- read_csv('Data/ICAR_NT.csv') %>%
    dplyr::select(`Participant Public ID`,
                  Response,
                  Attempt,
                  ANSWER, display, `Reaction Time`) %>%
    filter(Attempt != is.na(Attempt),
           display != 'instructions') %>%
    dplyr::select(-Attempt, -display) %>%
    mutate(Response = str_remove(Response, ".*_"),
           Response = str_remove(Response, ".png"),
           Correct = ifelse(Response == ANSWER, 1, 0)) %>%
    rename(ID = 1) %>%
    group_by(ID) %>%
    mutate(ICAR = sum(Correct)) %>%
    dplyr::select(ID, ICAR) %>%
  distinct()

nt_join_ICAR_mrat <- left_join(nt_onl %>%
            as.data.frame() %>%
            filter(str_detect(name, 'Mratio')) %>%
            slice(1:84) %>%
            mutate(ID = ID_nt_mrat$ID),
          ICAR_NT %>%
            filter(ID %in% ID_nt_mrat$ID),
          by = 'ID')

nt_join_ICAR_meta_d <- left_join(nt_onl %>%
            as.data.frame() %>%
            filter(str_detect(name, 'meta_d')) %>%
            slice(1:84) %>%
            mutate(ID = ID_nt_mrat$ID),
          ICAR_NT %>%
            filter(ID %in% ID_nt_mrat$ID),
          by = 'ID')

nt_join_ICAR_d1 <- left_join(nt_onl %>%
            as.data.frame() %>%
            filter(str_detect(name, 'd1')) %>%
            slice(1:84) %>%
            mutate(ID = ID_nt_mrat$ID),
          ICAR_NT %>%
            filter(ID %in% ID_nt_mrat$ID),
          by = 'ID')

cor.test(nt_join_ICAR_mrat$mean, nt_join_ICAR$ICAR)
cor.test(nt_join_ICAR_meta_d$mean, nt_join_ICAR$ICAR)
cor.test(nt_join_ICAR_d1$mean, nt_join_ICAR$ICAR)

# ── 14. SUPPLEMENT: EXCLUSION ANALYSIS -------------------------------

### Experiment 1 ------

run_lmer(correct ~ kmed + (1|ID), lmer_data_exc %>% filter(group == 'CCD (Online)'))
run_lmer(correct ~ kmed + (1|ID), lmer_data     %>% filter(group == 'NT (Online)'))
run_lmer(correct ~ kmed * group + (1|ID), lmer_data_exc %>% filter(group %in% c('NT (Online)', 'CCD (Online)')))

run_lmer(conf ~ kmed + (1|ID), lmer_data_exc %>% filter(group == 'CCD (Online)'))
run_lmer(conf ~ kmed + (1|ID), lmer_data     %>% filter(group == 'NT (Online)'))
run_lmer(conf ~ kmed * group + (1|ID), lmer_data_exc %>% filter(group %in% c('NT (Online)', 'CCD (Online)')))

### Experiment 3 ---------

run_lmer(Correct  ~ kmed * group + (1|ID), data_vr_tests %>% filter(!ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(conf  ~ kmed * group + (1|ID), data_vr_tests %>% filter(!ID %in% c('24104011', '24105011', '24109011', '24112011')))

#### BINO ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular'))

#### LAT ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized'))

#### MONO ----

run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(Correct  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular'))

run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
run_lmer(conf  ~ kmed + (1|ID), data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular'))

#### Computational -----------------------------------------------------------

outputONL_CCD_exc <- readRDS('HMetaDFit/Mratio/ONL_CCD_exc.rdata')

Results_ONL_CCD_exc <- HMeta_post_clean(outputONL_CCD_exc, 'CCD (ONL)', F)

data_for_permute_exp1_CCD_exc  <- checkBothRDK_w_exc %>% filter(group=='CCD (Online)') %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)
onl_CCDvsNTp_exc  <- permute_metad_group_mrat(data_for_permute_exp1_CCD_exc, data_for_permute_exp1_NT, cores = 10, nreps = nreps);

data_for_permute_exp3_CCD_exc  <- data_vr %>% filter(BinocularTrial   == 1, Group=='CCD', !PID %in% exc_ids)
data_for_permute_exp3_CCD_L_exc<- data_vr %>% filter(LateralizedTrial == 1, Group=='CCD', !PID %in% exc_ids)
data_for_permute_exp3_CCD_M_exc<- data_vr %>% filter(MonocularTrial   == 1, Group=='CCD', !PID %in% exc_ids)

vr_CCDvsNTp_exc  <- permute_metad_group_vr(data_for_permute_exp3_CCD_exc,   data_for_permute_exp3_NT,cores = 10, nreps = nreps);
vr_CCDvsNTp_L_exc<- permute_metad_group_vr(data_for_permute_exp3_CCD_L_exc, data_for_permute_exp3_NT_L,cores = 10, nreps = nreps);
vr_CCDvsNTp_M_exc<- permute_metad_group_vr(data_for_permute_exp3_CCD_M_exc, data_for_permute_exp3_NT_M,cores = 10, nreps = nreps);
