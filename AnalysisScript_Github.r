# RDK Script --------------------------------------------------------------

# Load libraries
library(ggpubr)
library(tidyverse)
library(patchwork)
library(doParallel)
library(rstanarm)
library(brms)
library(magrittr)
library(reshape2)
library(rjags)
library(coda)
library(lattice)
library(lme4)
library(broom)
library(ggmcmc)
library(foreach)
library(sjPlot)
library(webshot)
library(tidyplots)

source('RDKUtilityFuncs.R')

# Load Data ---------------------------------------------------------------

au_CCDRDKcheck <- read.csv('Data/Experiment1_CCD.csv')
online_NTRDKcheck <- read.csv('Data/Experiment1_NT.csv')
lab_NTRDKcheck <- read.csv('Data/Experiment2_CCD.csv')
us_CCDRDKcheck <- read.csv('Data/Experiment2_NT.csv')
data_vr <- read.csv('Data/Experiment3.csv')

# Plot diagnostics -------------------------------------------------------------
colpal <- c('#931F1D', '#FE938C', '#475B5A', '#81D6E3')

# Plot calibration coherence
(us_CCDRDKcheck %>% check_coherence(groupname = 'CCD (Online)', legend = F)+
au_CCDRDKcheck %>% check_coherence('CCD (MRI)', legend = F))/
(online_NTRDKcheck %>% check_coherence('NT (Online)', legend = F)+
lab_NTRDKcheck %>% check_coherence('NT (MRI)', legend = F))

au_CCDRDKcheck %>% filter(ID %in% c('1001', '1027')) %>% check_coherence('CCD (MRI)', legend = F)

# Check participant confidence reporting consistency
us_CCDRDKcheck %>% filter(conf != 0, type == 'main', mistake %in% c(0, NA, FALSE)) %>%
  ggplot(aes(conf, ID, fill = ID, group = ID)) +
   ggridges::geom_density_ridges(scale = 1, jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.2) +
   labs(x = 'Confidence', y = 'ID') +
   theme_bw() +
   theme(text = element_text(size = 18), legend.position = 'none')

au_CCDRDKcheck %>% filter(conf != 0, type == 'main', mistake %in% c(0, NA, FALSE)) %>%
  ggplot(aes(conf, ID, fill = ID, group = ID)) +
  ggridges::geom_density_ridges(scale = 1, jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.2) +
  labs(x = 'Confidence', y = 'ID') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')

online_NTRDKcheck %>% filter(conf != 0, type == 'main', mistake %in% c(0, NA, FALSE)) %>%
  ggplot(aes(conf, ID, fill = ID, group = ID)) +
  ggridges::geom_density_ridges(scale = 1, jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.2) +
  labs(x = 'Confidence', y = 'ID') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')

lab_NTRDKcheck %>% filter(conf != 0, type == 'main', mistake %in% c(0, NA, FALSE)) %>%
  ggplot(aes(conf, ID, fill = ID, group = ID)) +
  ggridges::geom_density_ridges(scale = 1, jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.2) +
  labs(x = 'Confidence', y = 'ID') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')

# Check the number of times participants reported they had made a mistake
us_CCDRDKcheck %>%
  filter(conf != 0, type == 'main') %>%
  group_by(ID) %>%
  mutate(mistake_sum = ifelse(mistake %in% c(T), 1, 0), mistake_sum = sum(mistake_sum)) %>%
  select(mistake_sum, ID) %>%
  unique() %>%
  ggplot(aes(ID, mistake_sum, fill = ID)) +
   geom_col()+
   coord_cartesian(ylim = c(1, 25))+
   labs(x = 'ID', y = 'Mistake') +
   theme_bw() +
   theme(text = element_text(size = 18), legend.position = 'none')

au_CCDRDKcheck %>%
  filter(conf != 0, type == 'main') %>%
  group_by(ID) %>%
  mutate(mistake_sum = ifelse(mistake %in% c(T), 1, 0), mistake_sum = sum(mistake_sum)) %>%
  select(mistake_sum, ID) %>%
  unique() %>%
  ggplot(aes(ID, mistake_sum, fill = ID)) +
  geom_col()+
  coord_cartesian(ylim = c(1, 25))+
  labs(x = 'ID', y = 'Mistake') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')

online_NTRDKcheck %>%
  filter(conf != 0, type == 'main') %>%
  group_by(ID) %>%
  mutate(mistake_sum = ifelse(mistake %in% c(T), 1, 0), mistake_sum = sum(mistake_sum)) %>%
  select(mistake_sum, ID) %>%
  unique() %>%
  ggplot(aes(ID, mistake_sum, fill = ID)) +
  geom_col()+
  coord_cartesian(ylim = c(1, 25))+
  labs(x = 'ID', y = 'Mistake') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')

lab_NTRDKcheck %>%
  filter(conf != 0, type == 'main') %>%
  group_by(ID) %>%
  mutate(mistake_sum = ifelse(mistake %in% c(T), 1, 0), mistake_sum = sum(mistake_sum)) %>%
  select(mistake_sum, ID) %>%
  unique() %>%
  ggplot(aes(ID, mistake_sum, fill = ID)) +
  coord_cartesian(ylim = c(1, 25))+
  geom_col()+
  labs(x = 'ID', y = 'Mistake') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = 'none')


# Remove participants if any calibration check failed --------------------------
# Also remove trials in which participant declared they made a mistake

# 2004: Exceeded coherence threshold of 0.5
# 5f20a850b6d5ed3d186bc3c3: Exceed coherence threshold of 0.5
# 1001: Exceed coherence threshold of 0.5
# 1004: Self-reported mistakes were not marked
# 1006: Reported confidence values had no variance and reported number of mistakes made was above 20
# 1027: Reported number of mistakes made was above 20
# 3002: Reported confidence values had no variance

ONL_CCDRDKcheck    <- us_CCDRDKcheck %>% filter(ID != '2004')
ONL_NTRDKcheck     <- online_NTRDKcheck %>% filter(ID != '5f20a850b6d5ed3d186bc3c3')
LAB_CCDRDKcheck    <- au_CCDRDKcheck %>% filter(ID != '1001', ID != '1006', ID != '1027')
LAB_NTRDKcheck     <- lab_NTRDKcheck %>% filter(ID != '3002')

ONL_CCDRDKcheck  <- ONL_CCDRDKcheck  %>% filter(mistake == FALSE)
ONL_NTRDKcheck   <- ONL_NTRDKcheck   %>% filter(mistake == FALSE)
LAB_CCDRDKcheck  <- LAB_CCDRDKcheck  %>% filter(mistake == FALSE)
LAB_NTRDKcheck   <- LAB_NTRDKcheck   %>% filter(mistake == FALSE)

# Final group counts
length(unique(us_CCDRDKcheck$ID))
length(unique(online_NTRDKcheck$ID))
length(unique(au_CCDRDKcheck$ID))
length(unique(lab_NTRDKcheck$ID))

unique(au_CCDRDKcheck$ID)
length(unique(au_CCDRDKcheck$ID))

# Combine data sets ------------------------------------------------------------
checkBothRDK <- rbind(ONL_CCDRDKcheck  %>% mutate(correct = as.numeric(correct)),
                      ONL_NTRDKcheck   %>% mutate(correct = as.numeric(correct)),
                      LAB_CCDRDKcheck  %>% mutate(correct = as.numeric(correct)),
                      LAB_NTRDKcheck   %>% mutate(correct = as.numeric(correct)))

checkBothRDK$highCoh <- ifelse(checkBothRDK$kmed == "kmed x 2", 1, 0)

# Determine if the testing environment potentially impacted performance --------
accuracy_data_vr <- data_vr %>%
  rename(ID=PID, environment=Environment) %>%
  filter(TrialType == 'Main',
         Presentation == 'Binocular') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Accuracy=mean(Correct))

accuracy_data_vr$ID <- as.character(accuracy_data_vr$ID)

accuracy_data_non_vr <- checkBothRDK %>%
  filter(type == 'main') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Accuracy=mean(correct))

accuracy_data <- rbind(accuracy_data_vr, accuracy_data_non_vr)

res = aov(Accuracy ~ environment + CCD, data=accuracy_data);
summary(res)

title <- "Effect of environment on participant accuracy"
model <- lm(formula = Accuracy ~ Environment + CCD, data=accuracy_data)
plot_model(model, title=title, show.values=TRUE, show.p=TRUE)
tab_model(model, show.re.var=TRUE, dv.labels=title, file="Stats/participant_accuracy_lmer.html")
webshot("Stats/participant_accuracy_lmer.html", "Stats/participant_accuracy_lmer.png")

data_vr$Confidence <- scale(data_vr$ReportedVeryConfident)
confidence_data_vr <- data_vr %>%
  filter(TrialType == 'Main') %>%
  group_by(PID, Environment, CCD) %>%
  summarise(Confidence=mean(Confidence))

confidence_data_vr$PID <- as.character(confidence_data_vr$PID)

checkBothRDKConfidenceTrialsOnly <- checkBothRDK %>% filter(conf != 0)
checkBothRDKConfidenceTrialsOnly$Confidence <- scale(checkBothRDKConfidenceTrialsOnly$conf)

confidence_data_non_vr <- checkBothRDKConfidenceTrialsOnly %>%
  filter(type == 'main') %>%
  group_by(ID, environment, CCD) %>%
  summarise(Confidence=mean(Confidence))

confidence_data_non_vr <- confidence_data_non_vr %>%
  rename(
    PID = ID,
    Environment = environment
  )

confidence_data <- rbind(confidence_data_vr, confidence_data_non_vr)

res = aov(Confidence ~ Environment + CCD, data=confidence_data);
summary(res)

# Calculate base k (VR) --------------------------------------------------------
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

# Calculate base k (non-VR) ----------------------------------------------------
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

# General group-wise plots for Exp1 and Exp 2-----------------------------------------------------

# Plot group reaction times vs high / low coherence and correct / incorrect
CCDcheck1 <- ggplot(checkBothRDK %>%
                     filter(type == 'main', kmed != 'kmed', RT < 5000, mistake %in% c(0, NA, FALSE)) %>%
                     mutate(correct = ifelse(correct == 1, 'Correct: Yes', 'Correct:No')))+
  geom_density(aes(RT, fill = group), show.legend = T, alpha = 0.7)+
  scale_fill_manual(values = colpal)+
  facet_wrap(kmed ~ correct)+
  labs(x = 'Reaction Time (ms)')+
  theme_minimal()+
  theme(axis.title.y = element_blank())
CCDcheck1

# Conduct non-parametric statistical testing of the group reaction time medians
kruskal.test(RT ~ group, data = rbind(ONL_CCDRDKcheck, ONL_NTRDKcheck))
kruskal.test(RT ~ group, data = rbind(LAB_CCDRDKcheck, LAB_NTRDKcheck))
kruskal.test(RT ~ group, data = rbind(LAB_CCDRDKcheck, ONL_CCDRDKcheck))
kruskal.test(RT ~ group, data = rbind(ONL_NTRDKcheck, LAB_NTRDKcheck))

# Calculate the mean number of correct responses for the two difference coherence levels
av_cor <- checkBothRDK %>%
  group_by(group, ID, kmed) %>%
  filter(type == 'main', kmed != 'kmed', mistake %in% c(0, NA, FALSE)) %>%
  mutate(av_cor = sum(corAdjusted)/length(corAdjusted),
         kmed = ifelse(kmed == 'kmed x 2', 'x2', 'x0.5')) %>%
  group_by(group, kmed) %>%
  dplyr::select(ID, av_cor, kmed) %>%
  distinct()

# Plot the mean number of correct responses for the two difference coherence levels
CCDcheck3 <- ggplot(av_cor, aes(kmed, av_cor, color = group, group = group))+
  geom_jitter(shape =21, width = 0.1, alpha = 0.1, aes(fill = group), colour = 'black')+
  stat_summary()+
  stat_summary(geom = 'line', show.legend = F)+
  # stat_compare_means(aes(kmed, corAdjusted, color = group, group = group),
                     # label = 'p.format', label.y = 0.55, method = 'kruskal.test', show.legend = F)+
  labs(x = 'Coherence', y = 'Correct'#, title = 'Correct Answers in Main'
       )+
  scale_color_manual(values = colpal)+
  scale_fill_manual(values = colpal)+
  coord_cartesian(ylim = c(0.5, 1))+
  theme_minimal()+
  theme(legend.position = 'none')
CCDcheck3

# Calculate the mean confidence rating for the two difference coherence levels
av_co <- checkBothRDK %>%
  group_by(group, ID, kmed) %>%
  filter(type == 'main', kmed != 'kmed', conf != 0, mistake %in% c(0, NA, FALSE)) %>%
  mutate(av_co = sum(conf)/length(conf),
         kmed = ifelse(kmed == 'kmed x 2', 'x2', 'x0.5')) %>%
  group_by(group, kmed) %>%
  dplyr::select(ID, av_co, kmed) %>%
  distinct()

# Plot the mean confidence rating for the two difference coherence levels
CCDcheck4<- ggplot(av_co, aes(kmed, av_co, color = group, group = group))+
  geom_jitter(shape =21, width = 0.1, alpha = 0.1, aes(fill = group), colour = 'black', show.legend = F)+
  stat_summary()+
  stat_summary(geom = 'line', show.legend = F)+
  labs(x = 'Coherence', y = 'Confidence'#,
       #title = 'Confidence by Coherence'
       )+
  scale_color_manual(values = colpal, name = 'Group')+
  scale_fill_manual(values = colpal)+
  coord_cartesian(ylim = c(50, 100))+
  theme_minimal()+
  theme(legend.position = 'none',
        legend.justification.top = 2,
        legend.background = element_rect(color = 'black'),
        legend.title = element_blank())
CCDcheck4

# Plot the mean accuracy and confidence rating for the two difference coherence levels
(CCDcheck3 | CCDcheck4)  &
  theme(text = element_text(size = 18),
        panel.grid = element_blank())

# Conduct non-parametric statistical testing of the group accuracy medians
t.test(av_cor ~ kmed, data = av_cor %>% filter(group == 'CCD (Online)'), paired = T) #exp1
t.test(av_cor ~ kmed, data = av_cor %>% filter(group == 'NT (Online)'), paired = T) #exp1
t.test(av_cor ~ kmed, data = av_cor %>% filter(group == 'CCD (MRI)'), paired = T) #exp2
t.test(av_cor ~ kmed, data = av_cor %>% filter(group == 'NT (MRI)'), paired = T) #exp2

aov(av_cor ~ group*kmed, data = av_cor %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
aov(av_cor ~ group*kmed, data = av_cor %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% summary()

# Conduct non-parametric statistical testing of the group confidence medians
t.test(av_co ~ kmed, data = av_co %>% filter(group == 'NT (Online)'), paired = T) #exp1
t.test(av_co ~ kmed, data = av_co %>% filter(group == 'CCD (Online)'), paired = T) #exp1
t.test(av_co ~ kmed, data = av_co %>% filter(group == 'NT (MRI)'), paired = T) #exp2
t.test(av_co ~ kmed, data = av_co %>% filter(group == 'CCD (MRI)'), paired = T) #exp2

aov(av_co ~ group*kmed, data = av_co %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
aov(av_co ~ group*kmed, data = av_co %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% summary()

## Exp 1 & 2 key stat tests ----

# Conduct LMER assessment
lmer_data <- checkBothRDK %>%
  group_by(group, ID, kmed) %>%
  filter(type == 'main',
         #kmed != 'kmed',
         mistake %in% c(0, NA, FALSE),
         conf != 0)

### Exp 1 ----

lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)')) %>% summary()
lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)')) %>% effectsize::effectsize()

lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)', !ID %in% c('2008', '2011'))) %>% summary()

lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (Online)')) %>% summary()
lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (Online)')) %>% effectsize::effectsize()

lmerTest::lmer(correct  ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
lmerTest::lmer(correct  ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)')) %>% summary()
lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)')) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (Online)', !ID %in% c('2008', '2011'))) %>% summary()

lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (Online)')) %>% summary()
lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (Online)')) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% summary()
lmerTest::lmer(conf ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (Online)', 'CCD (Online)'))) %>% effectsize::effectsize()

### Exp 2 ----

lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (MRI)')) %>% summary()
lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (MRI)')) %>% effectsize::effectsize()

lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (MRI)')) %>% summary()
lmerTest::lmer(correct  ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (MRI)')) %>% effectsize::effectsize()

lmerTest::lmer(correct  ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% summary()
lmerTest::lmer(correct  ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (MRI)')) %>% summary()
lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='CCD (MRI)')) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (MRI)')) %>% summary()
lmerTest::lmer(conf ~ kmed + (1|ID), data = lmer_data %>% filter(group=='NT (MRI)')) %>% effectsize::effectsize()

lmerTest::lmer(conf ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% summary()
lmerTest::lmer(conf ~ kmed * group + (1|ID), data = lmer_data %>% filter(group %in% c('NT (MRI)', 'CCD (MRI)'))) %>% effectsize::effectsize()

# General group-wise plots for Exp 3 -----------------------------------------------------

colpal_vr <- c('#B80C09', '#0E1C36')

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

# Calculate the mean number of correct responses for the two difference coherence levels
av_cor_vr <- data_vr_tests %>%
  group_by(group, ID, kmed, PresentationType) %>%
  filter(type == 'Main', kmed != 'kmed') %>%
  mutate(av_cor = sum(Correct)/length(Correct)) %>%
  group_by(group, kmed, PresentationType) %>%
  dplyr::select(ID, group, av_cor, kmed, PresentationType) %>%
  distinct()

# Plot the mean number of correct responses for the two difference coherence levels
vr_cor_plot <- ggplot(av_cor_vr, aes(kmed, av_cor, color = group, group = group))+
  geom_jitter(shape =21, width = 0.1, alpha = 0.1, aes(fill = group), colour = 'black')+
  stat_summary()+
  stat_summary(geom = 'line', show.legend = F)+

  labs(x = 'Coherence', y = 'Correct'
       )+
  facet_wrap(~PresentationType)+
  scale_color_manual(values = colpal_vr)+
  scale_fill_manual(values = colpal_vr)+
  coord_cartesian(ylim = c(0.5, 1))+
  theme_minimal()+
  theme(legend.position = 'none',
        axis.title.x = element_blank())
vr_cor_plot

# Calculate the mean confidence rating for the two difference coherence levels
av_co_vr <- data_vr_tests %>%
  group_by(group, ID, kmed, PresentationType) %>%
  filter(type == 'Main', kmed != 'kmed') %>%
  mutate(av_co = sum(conf)/length(conf)) %>%
  group_by(group, kmed, PresentationType) %>%
  dplyr::select(ID, group, av_co, kmed, PresentationType) %>%
  distinct()

# Plot the mean confidence rating for the two difference coherence levels
vr_co_plot <- ggplot(av_co_vr, aes(kmed, av_co, color = group, group = group))+
  geom_jitter(shape =21, width = 0.1, alpha = 0.1, aes(fill = group), colour = 'black', show.legend = F)+
  stat_summary()+
  stat_summary(geom = 'line', show.legend = F)+
  labs(x = 'Coherence', y = 'Confidence'#,
       #title = 'Confidence by Coherence'
       )+
  facet_wrap(~PresentationType)+
  scale_color_manual(values = colpal_vr, name = 'Group')+
  scale_fill_manual(values = colpal_vr)+
  theme_minimal()+
  theme(legend.position = 'none',
        legend.justification.top = 2,
        legend.background = element_rect(color = 'black'),
        legend.title = element_blank(),
        strip.text.x = element_blank())
vr_co_plot

# Plot the mean accuracy and confidence rating for the two difference coherence levels
(vr_cor_plot / vr_co_plot)  &
  theme(text = element_text(size = 18),
        panel.grid = element_blank())

# Conduct non-parametric statistical testing of the group accuracy medians
t.test(av_cor ~ kmed, data = av_cor_vr %>% filter(group == 'CCD', PresentationType=='Binocular'))
t.test(av_cor ~ kmed, data = av_cor_vr %>% filter(group == 'CCD', PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
t.test(av_cor ~ kmed, data = av_cor_vr %>% filter(group == 'NT',  PresentationType=='Binocular'))
aov(av_cor ~ kmed*group, data = av_cor_vr %>% filter(PresentationType=='Binocular')) %>% summary
aov(av_cor ~ kmed*group, data = av_cor_vr %>% filter(PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011'))) %>% summary

# Conduct non-parametric statistical testing of the group confidence medians
t.test(av_co ~ kmed, data = av_co_vr %>% filter(group == 'CCD', PresentationType=='Binocular'))
t.test(av_co ~ kmed, data = av_co_vr %>% filter(group == 'CCD', PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011')))
t.test(av_co ~ kmed, data = av_co_vr %>% filter(group == 'NT',  PresentationType=='Binocular'))
aov(av_co ~ kmed*group, data = av_co_vr %>% filter(PresentationType=='Binocular')) %>% summary
aov(av_co ~ kmed*group, data = av_co_vr %>% filter(PresentationType=='Binocular', !ID %in% c('24104011', '24105011', '24109011', '24112011'))) %>% summary

## Exp 3 key stat tests ----

lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests)  %>% summary()
lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests)  %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(!ID %in% c('24104011', '24105011', '24109011', '24112011')))  %>% summary()

lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests)  %>% summary()
lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests)  %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(!ID %in% c('24104011', '24105011', '24109011', '24112011')))  %>% summary()

### bino ----

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Binocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Binocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Binocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Binocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Binocular')) %>% effectsize::effectsize()

### lat ----

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Lateralized')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Lateralized')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Lateralized')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Lateralized')) %>% summary()
lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Lateralized')) %>% effectsize::effectsize()

### mono ----

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular')) %>% effectsize::effectsize()

lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(Correct  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Monocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='CCD', PresentationType=='Monocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed + (1|ID), data = data_vr_tests %>% filter(group=='NT', PresentationType=='Monocular')) %>% effectsize::effectsize()

lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Monocular')) %>% summary()
lmerTest::lmer(conf  ~ kmed * group + (1|ID), data = data_vr_tests %>% filter(PresentationType=='Monocular')) %>% effectsize::effectsize()

## Reaction times -----

summary(lm(ReactionTime ~ PresentationType * Group, data_vr))

# H-Meta D Ratio calculation preparation ----------------------------------------------------------

# Create data set for LAB CCD participants
LAB_CCD_metaDat <- HMetaGroupPrep(checkBothRDK, 'CCD (MRI)')

# Create data set for ONL CCD participants
ONL_CCD_metaDat <- HMetaGroupPrep(checkBothRDK, 'CCD (Online)')

# Create data set for VR CCD participants
VR_CCD_metaDat  <- HMetaGroupPrep_vr(data_vr %>% filter(BinocularTrial == 1), 'CCD')
VR_CCD_metaDat_L<- HMetaGroupPrep_vr(data_vr %>% filter(LateralizedTrial == 1), 'CCD')
VR_CCD_metaDat_M<- HMetaGroupPrep_vr(data_vr %>% filter(MonocularTrial == 1), 'CCD')

# Create combined data set for the CCD participants
all_CCD_metaDat1 <- cbind(LAB_CCD_metaDat[[1]], ONL_CCD_metaDat[[1]]); colnames(all_CCD_metaDat1) <- 1:21
all_CCD_metaDat2 <- cbind(LAB_CCD_metaDat[[2]], ONL_CCD_metaDat[[2]]); colnames(all_CCD_metaDat2) <- 1:21

all_CCD_metaDat  <- list(S1mat = all_CCD_metaDat1, S2mat = all_CCD_metaDat2)

# Create data set for AU NT participants
ONL_NT_metaDat <- HMetaGroupPrep(checkBothRDK, 'NT (Online)')

# Create data set for  CCD participants
LAB_NT_metaDat <- HMetaGroupPrep(checkBothRDK, 'NT (MRI)')

# Create data set for VR CCD participants
VR_NT_metaDat  <- HMetaGroupPrep_vr(data_vr %>% filter(BinocularTrial == 1), 'NT')
VR_NT_metaDat_L<- HMetaGroupPrep_vr(data_vr %>% filter(LateralizedTrial == 1), 'NT')
VR_NT_metaDat_M<- HMetaGroupPrep_vr(data_vr %>% filter(MonocularTrial == 1), 'NT')

# Hierarchical meta_d group function ------------------------------------------------------

setwd('HMeta-d-master/R/')
source('Function_metad_group.R')

# Fit all data
#outputAU_CCD     <- metad_group(nR_S1 = list(LAB_CCD_metaDat[[1]]), nR_S2 = list(LAB_CCD_metaDat[[2]]));
#outputUS_CCD     <- metad_group(nR_S1 = list(ONL_CCD_metaDat[[1]]), nR_S2 = list(ONL_CCD_metaDat[[2]]));
#outputVR_CCD     <- metad_group(nR_S1 = list(VR_CCD_metaDat[[1]]),  nR_S2 = list(VR_CCD_metaDat[[2]]));
#outputVR_CCD_L   <- metad_group(nR_S1 = list(VR_CCD_metaDat_L[[1]]),nR_S2 = list(VR_CCD_metaDat_L[[2]]));
#outputVR_CCD_M   <- metad_group(nR_S1 = list(VR_CCD_metaDat_M[[1]]),nR_S2 = list(VR_CCD_metaDat_M[[2]]));
#outputOnline_NT  <- metad_group(nR_S1 = list(ONL_NT_metaDat[[1]]),  nR_S2 = list(ONL_NT_metaDat[[2]]));
#outputLab_NT     <- metad_group(nR_S1 = list(LAB_NT_metaDat[[1]]),  nR_S2 = list(LAB_NT_metaDat[[2]]));
#outputVR_NT      <- metad_group(nR_S1 = list(VR_NT_metaDat[[1]]),   nR_S2 = list(VR_NT_metaDat[[2]]));
#outputVR_NT_L    <- metad_group(nR_S1 = list(VR_NT_metaDat_L[[1]]), nR_S2 = list(VR_NT_metaDat_L[[2]]));
#outputVR_NT_M    <- metad_group(nR_S1 = list(VR_NT_metaDat_M[[1]]), nR_S2 = list(VR_NT_metaDat_M[[2]]));
#
## Save the output
#saveRDS(outputAU_CCD,    'HMetaDFit/LAB_CCD.rdata')
#saveRDS(outputUS_CCD,    'HMetaDFit/ONL_CCD.rdata')
#saveRDS(outputVR_CCD,    'HMetaDFit/VR_CCD.rdata')
#saveRDS(outputVR_CCD_L,  'HMetaDFit/VR_CCD_L.rdata')
#saveRDS(outputVR_CCD_M,  'HMetaDFit/VR_CCD_M.rdata')
#saveRDS(outputLab_NT,    'HMetaDFit/LAB_NT.rdata')
#saveRDS(outputOnline_NT, 'HMetaDFit/ONL_NT.rdata')
#saveRDS(outputVR_NT,     'HMetaDFit/VR_NT.rdata')
#saveRDS(outputVR_NT_L,   'HMetaDFit/VR_NT_L.rdata')
#saveRDS(outputVR_NT_M,   'HMetaDFit/VR_NT_M.rdata')

outputLAB_CCD   <- readRDS('HMetaDFit/LAB_CCD.rdata')
outputONL_CCD   <- readRDS('HMetaDFit/ONL_CCD.rdata')
outputVR_CCD    <- readRDS('HMetaDFit/VR_CCD.rdata')
outputVR_CCD_L  <- readRDS('HMetaDFit/VR_CCD_L.rdata')
outputVR_CCD_M  <- readRDS('HMetaDFit/VR_CCD_M.rdata')
outputLAB_NT    <- readRDS('HMetaDFit/LAB_NT.rdata')
outputONL_NT    <- readRDS('HMetaDFit/ONL_NT.rdata')
outputVR_NT     <- readRDS('HMetaDFit/VR_NT.rdata')
outputVR_NT_L   <- readRDS('HMetaDFit/VR_NT_L.rdata')
outputVR_NT_M   <- readRDS('HMetaDFit/VR_NT_M.rdata')

# Extract the cleaned up values
Results_LAB_CCD <- HMeta_post_clean(outputLAB_CCD, 'CCD (LAB)', F)
Results_ONL_CCD <- HMeta_post_clean(outputONL_CCD, 'CCD (ONL)', F)
Results_VR_CCD  <- HMeta_post_clean(outputVR_CCD, 'CCD (VR)', F)
Results_VR_CCD_L<- HMeta_post_clean(outputVR_CCD_L, 'CCD (VR L)', F)
Results_VR_CCD_M<- HMeta_post_clean(outputVR_CCD_M, 'CCD (VR M)', F)
Results_ONL_NT  <- HMeta_post_clean(outputONL_NT, 'NT (ONL)', F)
Results_LAB_NT  <- HMeta_post_clean(outputLAB_NT, 'NT (LAB)', F)
Results_VR_NT   <- HMeta_post_clean(outputVR_NT, 'NT (VR)', F)
Results_VR_NT_L <- HMeta_post_clean(outputVR_NT_L, 'NT (VR L)', F)
Results_VR_NT_M <- HMeta_post_clean(outputVR_NT_M, 'NT (VR M)', F)

# Extract the Markov Chain Monte Carlo (MCMC) values
mcmcGroup       <- rbind(Results_LAB_CCD[[1]], Results_ONL_CCD[[1]],
                         Results_ONL_NT[[1]],  Results_LAB_NT[[1]])
mcmcGroup_vr    <- rbind(Results_VR_CCD[[1]], Results_VR_NT[[1]],
                         Results_VR_CCD_L[[1]], Results_VR_NT_L[[1]],
                         Results_VR_CCD_M[[1]], Results_VR_NT_M[[1]])
mcmcGroup2      <- rbind(Results_All_CCD[[1]], Results_All_NT[[1]])

# Generate posterior distributions for the mu M-ratio values
statGroup <- data.frame(
  CCD_LAB     = c(exp(Results_LAB_CCD[[3]]$mean[Results_LAB_CCD[[3]]$name == "mu_logMratio"]),   exp(Results_LAB_CCD[[3]]$lower[Results_LAB_CCD[[3]]$name == "mu_logMratio"]),    exp(Results_LAB_CCD[[3]]$upper[Results_LAB_CCD[[3]]$name == "mu_logMratio"]),
                  exp(Results_LAB_CCD[[3]]$mean[Results_LAB_CCD[[3]]$name == "sigma_logMratio"]),exp(Results_LAB_CCD[[3]]$lower[Results_LAB_CCD[[3]]$name == "sigma_logMratio"]), exp(Results_LAB_CCD[[3]]$upper[Results_LAB_CCD[[3]]$name == "sigma_logMratio"])),
  CCD_ONL     = c(exp(Results_ONL_CCD[[3]]$mean[Results_ONL_CCD[[3]]$name == "mu_logMratio"]),   exp(Results_ONL_CCD[[3]]$lower[Results_ONL_CCD[[3]]$name == "mu_logMratio"]),    exp(Results_ONL_CCD[[3]]$upper[Results_ONL_CCD[[3]]$name == "mu_logMratio"]),
                  exp(Results_ONL_CCD[[3]]$mean[Results_ONL_CCD[[3]]$name == "sigma_logMratio"]),exp(Results_ONL_CCD[[3]]$lower[Results_ONL_CCD[[3]]$name == "sigma_logMratio"]), exp(Results_ONL_CCD[[3]]$upper[Results_ONL_CCD[[3]]$name == "sigma_logMratio"])),
  NT_ONL      = c(exp(Results_ONL_NT[[3]]$mean[Results_ONL_NT[[3]]$name == "mu_logMratio"]),   exp(Results_ONL_NT[[3]]$lower[Results_ONL_NT[[3]]$name == "mu_logMratio"]),    exp(Results_ONL_NT[[3]]$upper[Results_ONL_NT[[3]]$name == "mu_logMratio"]),
                  exp(Results_ONL_NT[[3]]$mean[Results_ONL_NT[[3]]$name == "sigma_logMratio"]),exp(Results_ONL_NT[[3]]$lower[Results_ONL_NT[[3]]$name == "sigma_logMratio"]), exp(Results_ONL_NT[[3]]$upper[Results_ONL_NT[[3]]$name == "sigma_logMratio"])),
  NT_LAB      = c(exp(Results_LAB_NT[[3]]$mean[Results_LAB_NT[[3]]$name == "mu_logMratio"]),   exp(Results_LAB_NT[[3]]$lower[Results_LAB_NT[[3]]$name == "mu_logMratio"]),    exp(Results_LAB_NT[[3]]$upper[Results_LAB_NT[[3]]$name == "mu_logMratio"]),
                  exp(Results_LAB_NT[[3]]$mean[Results_LAB_NT[[3]]$name == "sigma_logMratio"]),exp(Results_LAB_NT[[3]]$lower[Results_LAB_NT[[3]]$name == "sigma_logMratio"]), exp(Results_LAB_NT[[3]]$upper[Results_LAB_NT[[3]]$name == "sigma_logMratio"])),
  CCD_All     = c(exp(Results_All_CCD[[3]]$mean[Results_All_CCD[[3]]$name == "mu_logMratio"]),   exp(Results_All_CCD[[3]]$lower[Results_All_CCD[[3]]$name == "mu_logMratio"]),    exp(Results_All_CCD[[3]]$upper[Results_All_CCD[[3]]$name == "mu_logMratio"]),
                  exp(Results_All_CCD[[3]]$mean[Results_All_CCD[[3]]$name == "sigma_logMratio"]),exp(Results_All_CCD[[3]]$lower[Results_All_CCD[[3]]$name == "sigma_logMratio"]), exp(Results_All_CCD[[3]]$upper[Results_All_CCD[[3]]$name == "sigma_logMratio"])),
  NT_All      = c(exp(Results_All_NT[[3]]$mean[Results_All_NT[[3]]$name == "mu_logMratio"]),   exp(Results_All_NT[[3]]$lower[Results_All_NT[[3]]$name == "mu_logMratio"]),    exp(Results_All_NT[[3]]$upper[Results_All_NT[[3]]$name == "mu_logMratio"]),
                  exp(Results_All_NT[[3]]$mean[Results_All_NT[[3]]$name == "sigma_logMratio"]),exp(Results_All_NT[[3]]$lower[Results_All_NT[[3]]$name == "sigma_logMratio"]), exp(Results_All_NT[[3]]$upper[Results_All_NT[[3]]$name == "sigma_logMratio"])),

  Variable= c('mu', 'lower mu', 'upper mu', 'sigma', 'lower sigma', 'upper sigma')
)

statGroup[,1:2] <- sapply(statGroup[,1:2], as.numeric)
statGroup

exp(Results_VR_CCD[[3]]$mean[Results_VR_CCD[[3]]$name == "mu_logMratio"])
exp(Results_VR_NT[[3]]$mean[Results_VR_NT[[3]]$name == "mu_logMratio"])

exp(Results_VR_CCD_L[[3]]$mean[Results_VR_CCD_L[[3]]$name == "mu_logMratio"])
exp(Results_VR_NT_L[[3]]$mean[Results_VR_NT_L[[3]]$name == "mu_logMratio"])

exp(Results_VR_CCD_M[[3]]$mean[Results_VR_CCD_M[[3]]$name == "mu_logMratio"])
exp(Results_VR_NT_M[[3]]$mean[Results_VR_NT_M[[3]]$name == "mu_logMratio"])

# Test the output
mcmcMu <- mcmcGroup %>%
  filter(Parameter == "mu_logMratio")
mcmcSigma <- mcmcGroup %>%
  filter(Parameter == "sigma_logMratio")

# Plot the mu distributions for the in-lab participants
p1 <- mcmcGroup %>%
  filter(Parameter == "mu_logMratio") %>%
  ggplot(aes(exp(value))) +
  geom_histogram(aes(fill = group), binwidth = 0.03, alpha = 0.9, position = 'dodge') +
  scale_fill_manual(values = colpal)+
  labs(y="Sample Count", x = expression(paste(mu["ratio"]))) +
  coord_cartesian(xlim = c(0, 1.5))+
  scale_y_continuous(expand = c(0,0))+
  ggdist::theme_tidybayes()+
  theme(axis.text.x = element_text(size = 20, color = 'black'),
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.70),
        legend.background = element_blank(),
        panel.background = element_rect(fill = 'white'))

#ggsave(p1, filename = "LindaFigure1.png", bg = "transparent")
p1

# Plot the mu distributions for the in-lab participants
p1b <- mcmcGroup_vr %>%
  filter(Parameter == "mu_logMratio",
         group %in% c('CCD (VR)', 'NT (VR)')) %>%
  ggplot(aes(exp(value))) +
  geom_histogram(aes(fill = group), binwidth = 0.03, alpha = 0.9, position = 'dodge') +
  scale_fill_manual(values = colpal_vr)+
  labs(y="Sample Count", x = expression(paste(mu["ratio"]))) +
  coord_cartesian(xlim = c(0, 1.5))+
  scale_y_continuous(expand = c(0,0))+
  ggdist::theme_tidybayes()+
  theme(axis.text = element_text(size = 20, color = 'black'),
        axis.title = element_text(size = 20, color = 'black'),
        legend.text = element_text(size = 16, color = 'black'),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(color = 'black'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'))

#ggsave(p1, filename = "LindaFigure1.png", bg = "transparent")
p1b

# Plot the mu distributions for the in-lab participants
p1c <- mcmcGroup_vr %>%
  filter(Parameter == "mu_logMratio",
         group %in% c('CCD (VR L)', 'NT (VR L)')) %>%
  ggplot(aes(exp(value))) +
  geom_histogram(aes(fill = group), binwidth = 0.03, alpha = 0.9, position = 'dodge') +
  scale_fill_manual(values = colpal_vr[1:2])+
  labs(y="Sample Count", x = expression(paste(mu["ratio"]))) +
  coord_cartesian(xlim = c(0, 1.5))+
  scale_y_continuous(expand = c(0,0))+
  ggdist::theme_tidybayes()+
  theme(axis.text = element_text(size = 20, color = 'black'),
        axis.title = element_text(size = 20, color = 'black'),
        legend.text = element_text(size = 16, color = 'black'),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(color = 'black'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'))

#ggsave(p1, filename = "LindaFigure1.png", bg = "transparent")
p1c

# Plot the mu distributions for the in-lab participants
p1d <- mcmcGroup_vr %>%
  filter(Parameter == "mu_logMratio",
         group %in% c('CCD (VR M)', 'NT (VR M)')) %>%
  ggplot(aes(exp(value))) +
  geom_histogram(aes(fill = group), binwidth = 0.03, alpha = 0.9, position = 'dodge') +
  scale_fill_manual(values = colpal_vr[1:2])+
  labs(y="Sample Count", x = expression(paste(mu["ratio"]))) +
  coord_cartesian(xlim = c(0, 1.5))+
  scale_y_continuous(expand = c(0,0))+
  ggdist::theme_tidybayes()+
  theme(axis.text = element_text(size = 20, color = 'black'),
        axis.title = element_text(size = 20, color = 'black'),
        legend.text = element_text(size = 16, color = 'black'),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(color = 'black'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'))

#ggsave(p1, filename = "LindaFigure1.png", bg = "transparent")
p1d

# m-ratio spread indiv ----------------------------------------------------

make_block <- function(group, type, v) {
  data.frame(
    group = group,
    type  = type,
    value = v,
    stringsAsFactors = FALSE
  )
}

dat_long <- rbind(
    do.call(rbind, list(
       make_block("CCD", "BIN", ccd_b[1:10, 2]),
       make_block("CCD", "MON", ccd_m[1:10, 2]),
       make_block("CCD", "LAT", ccd_l[1:10, 2]),

       make_block("NT", "BIN", nt_b[1:13, 2]),
       make_block("NT", "MON", nt_m[1:13, 2]),
       make_block("NT", "LAT", nt_l[1:13, 2])
     )),
   do.call(rbind, list(
            make_block("CCD", "BIN", ccd_b[11:20, 2]),
            make_block("CCD", "MON", ccd_m[11:20, 2]),
            make_block("CCD", "LAT", ccd_l[11:20, 2]),

            make_block("NT", "BIN", nt_b[14:26, 2]),
            make_block("NT", "MON", nt_m[14:26, 2]),
            make_block("NT", "LAT", nt_l[14:26, 2])
          )),
    do.call(rbind, list(
        make_block("CCD", "BIN", ccd_b[21:30, 2]),
        make_block("CCD", "MON", ccd_m[21:30, 2]),
        make_block("CCD", "LAT", ccd_l[21:30, 2]),

        make_block("NT", "BIN", nt_b[27:39, 2]),
        make_block("NT", "MON", nt_m[27:39, 2]),
        make_block("NT", "LAT", nt_l[27:39, 2])
      ))
) %>%
  mutate(ID = rep(c(rep(1:10, 3), rep(1:13, 3)), 3),
         parm = c(rep('mratio', 10*3+13*3), rep('d1', 10*3+13*3), rep('metad', 10*3+13*3)))

dat_long

ggplot(dat_long %>% filter(type == 'LAT', parm == 'mratio'),
       aes(y = group, x = value, fill = group))+
  geom_vline(xintercept = 1, size = 2, colour = 'grey')+
  geom_jitter(shape = 21, height = 0.1, size = 3, alpha = 0.8)+
  scale_fill_manual(values = colpal_vr[1:2])+
  facet_wrap(~type)+
  labs(x = expression(paste('Indiv. ', mu, '-ratio')))+
  stat_compare_means(paired = F, method = 't.test')+
  coord_cartesian(xlim = c(0, 1.5))+
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5), limits = c(0, 2.5), expand = c(0,0))+
  theme_bw(base_size = 24)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text = element_blank())

## Correlations of outcomes ------

library(dplyr)
library(tidyr)
library(ggplot2)

results <- list()

pairs <- list(
  c("BIN", "LAT"),
  c("BIN", "MON"),
  c("LAT", "MON")
)

groups_present <- sort(unique(dat_long$group))
scopes <- c("ALL", groups_present)

for (p in unique(dat_long$parm)) {
  for (g in scopes) {

    # Filter
    df <- dat_long %>%
      filter(parm == p) %>%
      { if (g != "ALL") filter(., group == g) else . }

    # Align values by index within type
    df_aligned <- df %>%
      group_by(type) %>%
      mutate(idx = row_number()) %>%
      ungroup()

    wide <- df_aligned %>%
      select(parm, group, type, idx, value) %>%
      pivot_wider(names_from = type, values_from = value)

    for (pair in pairs) {
      a <- pair[1]; b <- pair[2]
      sub <- wide %>% select(all_of(c(a, b))) %>% drop_na()

      if (nrow(sub) > 1) {
        ct <- suppressWarnings(cor.test(sub[[a]], sub[[b]], method = "spearman", exact = FALSE))
        results[[length(results) + 1]] <- data.frame(
          parm       = p,
          group      = g,
          comparison = paste(a, b, sep = "-"),
          rho        = unname(ct$estimate),
          pval       = ct$p.value
        )
      } else {
        results[[length(results) + 1]] <- data.frame(
          parm       = p,
          group      = g,
          comparison = paste(a, b, sep = "-"),
          rho        = NA_real_,
          pval       = NA_real_
        )
      }
    }
  }
}

results_df <- bind_rows(results)

results_df <- results_df %>%
  mutate(
    sig = !is.na(pval) & pval <= 0.05,
    rho_lab = ifelse(is.na(rho), "", sprintf("rho=%.2f", rho)),
    group = factor(group, levels = c("ALL", sort(unique(as.character(group[group != "ALL"])))))
  )

ggplot() +
  geom_tile(
    data = subset(results_df, !sig | is.na(sig)),
    aes(x = group, y = comparison),
    fill = "white", color = "grey85"
  ) +
  geom_tile(
    data = subset(results_df, sig),
    aes(x = group, y = comparison, fill = rho),
    color = "white"
  ) +
  geom_text(
    data = results_df,
    aes(x = group, y = comparison, label = rho_lab),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limits = c(-1, 1), na.value = "white",
    name = "Spearman rho"
  ) +
  facet_wrap(~ parm, nrow = 1) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = 'none',
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

# d', meta-d', m-ratio spread -----------------------------------------------------------------

#LAB_CCD_metaDat
#ONL_CCD_metaDat
#VR_CCD_metaDat
do=0
if(do==1){
 ccd_onl     <- make_results(ONL_CCD_metaDat, 'CCD', 'ONL')
 ccd_lab     <- make_results(LAB_CCD_metaDat, 'CCD', 'LAB')
 ccd_b       <- make_results(VR_CCD_metaDat, "CCD", 'BINO')
 ccd_l       <- make_results(VR_CCD_metaDat_L, "CCD", 'LAT')
 ccd_m       <- make_results(VR_CCD_metaDat_M, "CCD", 'MONO')

 nt_onl      <- make_results(ONL_NT_metaDat, 'NT', 'ONL')
 nt_lab      <- make_results(LAB_NT_metaDat, 'NT', 'LAB')
 nt_b        <- make_results(VR_NT_metaDat, "NT", 'BINO')
 nt_l        <- make_results(VR_NT_metaDat_L, "NT", 'LAT')
 nt_m        <- make_results(VR_NT_metaDat_M, "NT", 'MONO')

meta_d_df <- rbind(ccd_onl, ccd_lab, ccd_b, ccd_l, ccd_m,
                   nt_onl, nt_lab, nt_b, nt_l, nt_m)

}

## PLOT LAB ONL --------

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'Mratio'),
                     !name %in% c('mu_logMratio', 'sigma_logMratio'),
                     type == 'ONL')
))

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'meta_d'),
                     !name %in% c('mu_meta_d'),
                     type == 'ONL')
))

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'd1'),
                     !name %in% c('mu_d1'),
                     type == 'ONL')
))

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'Mratio'),
                     !name %in% c('mu_logMratio', 'sigma_logMratio'),
                     type == 'LAB')
))

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'meta_d'),
                     !name %in% c('mu_meta_d'),
                     type == 'LAB')
))

summary(lm(mean ~ group,
           data =
            meta_d_df %>%
              filter(str_detect(name, 'd1'),
                     !name %in% c('mu_d1'),
                     type == 'LAB')
))

ggplot(meta_d_df %>%
         filter(str_detect(name, "mu_logMratio"),
                type %in% c('ONL', 'LAB')) %>%
         mutate(meld = paste(type, group)),
       aes(fct_rev(type), exp(mean), fill = meld))+
      #geom_jitter(data = meta_d_df %>%
      #              filter(str_detect(name, "d1"), type %in% c('ONL', 'LAB')) %>%
      #              mutate(meld = paste(type, group)),
      #            aes(x = type, y = abs(mean), fill = meld),
      #            shape = 21, alpha = 0.5, size = 2,
      #            position = position_dodge(width = 0.6))+
      stat_summary(data = meta_d_df %>%
                    filter(str_detect(name, "Mratio"), type %in% c('ONL', 'LAB')) %>%
                    mutate(meld = paste(type, group)),
                    aes(x = fct_rev(type), y = abs(mean)),
                    fun.data = "mean_cl_boot", linewidth = 2, size = 1, shape = NA,
                    position = position_dodge(width = 0.6))+
      stat_summary(shape=21, size = 1.2,
                  position = position_dodge(width = 0.6))+
      coord_cartesian(ylim=c(0, 3))+
      scale_y_continuous(expand = c(0,0), breaks = c(0, 1, 2, 3, 4, 5))+
      scale_fill_manual(values = colpal[c(1,3,2,4)])+
      theme_bw(base_size=18)&
      theme(legend.position = 'none',
            legend.title = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank())

meta_d_plot_1 <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_meta_d"),
                                           type %in% c('ONL', 'LAB')),
       aes(type, abs(mean), fill = group))

d_prime_plot_1 <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_d1"),
                                           type %in% c('ONL', 'LAB')),
       aes(type, abs(mean), fill = group))

m_ratio_plot_1 <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_logMratio"),
                                           type %in% c('ONL', 'LAB')),
       aes(type, exp(mean), fill = group))

(d_prime_plot_1|meta_d_plot_1|m_ratio_plot_1) &

  stat_summary(shape=21, size = 1.2)&
  coord_cartesian(ylim=c(0, 2.5))&
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1, 1.5, 2, 2.5))&
  scale_fill_manual(values = colpal_vr[1:2])&
  theme_bw(base_size=18)&
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

meta_d_df %>%
  mutate(name = str_remove(name, "\\s*\\[\\d+\\]$")) %>%
  filter(name %in% c('Mratio', 'd1', 'meta_d'),
         type == 'ONL') %>%

  ggplot(aes(abs(mean), fct_rev(name), fill = group))+
  geom_vline(xintercept = 1, size = 2, colour = 'grey')+
  geom_boxplot(alpha = 0.9, width = 0.6, outliers = F)+
  coord_cartesian(xlim = c(0, 2.5))+
  scale_x_continuous(breaks = c(0, 1, 2))+
  scale_fill_manual(values = colpal_vr[1:2])+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

## PLOT VR --------

meta_d_plot <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_meta_d"),
                                           type %in% c('BINO', 'MONO', 'LAT')),
       aes(type, abs(mean), fill = group))

d_prime_plot <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_d1"),
                                           type %in% c('BINO', 'MONO', 'LAT')),
       aes(type, abs(mean), fill = group))

m_ratio_plot <- ggplot(meta_d_df %>% filter(str_detect(name, "mu_logMratio"),
                                           type %in% c('BINO', 'MONO', 'LAT')),
       aes(type, exp(mean), fill = group))

patterns <- c("Mratio", "d1", "meta_d")


(d_prime_plot|meta_d_plot|m_ratio_plot) &

  stat_summary(geom = 'line', aes(group=group))&
  stat_summary(shape=21, size = 1.2)&
  coord_cartesian(ylim=c(0, 2))&
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1, 1.5, 2))&
  stat_compare_means(label = 'p.signif', label.y.npc = 0.3)&
  scale_fill_manual(values = colpal_vr[1:2])&
  theme_bw(base_size=18)&
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

meta_d_df %>%
  mutate(name = str_remove(name, "\\s*\\[\\d+\\]$")) %>%
  filter(name %in% c('Mratio', 'd1', 'meta_d'),
         type == 'LAT') %>%

  ggplot(aes(abs(mean), fct_rev(name), fill = group))+
  geom_vline(xintercept = 1, size = 2, colour = 'grey')+
  #geom_jitter(height = 0.1, shape = 21, size = 2, alpha = 0.2)+
  geom_boxplot(alpha = 0.9, width = 0.6, outliers = F)+
  coord_cartesian(xlim = c(0, 2.5))+
  scale_x_continuous(breaks = c(0, 1, 2))+
  scale_fill_manual(values = colpal_vr[1:2])+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

## Check by coherence ------

LAB_CCD_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'CCD (MRI)')
LAB_CCD_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'CCD (MRI)')

# Create data set for ONL CCD participants
ONL_CCD_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'CCD (Online)')
ONL_CCD_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'CCD (Online)')

# Create data set for AU NT participants
ONL_NT_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'NT (Online)')
ONL_NT_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'NT (Online)')

# Create data set for  CCD participants
LAB_NT_metaDat_hc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 2'),   'NT (MRI)')
LAB_NT_metaDat_lc <- HMetaGroupPrep(checkBothRDK %>% filter(kmed == 'kmed x 0.5'), 'NT (MRI)')

do=0
if(do==1){

 ccd_lab_hc <- make_results(LAB_CCD_metaDat_hc, 'CCD', 'LAB') %>% mutate(coh = 'x2')
 ccd_lab_lc <- make_results(LAB_CCD_metaDat_lc, 'CCD', 'LAB') %>% mutate(coh = 'x0.5')
 ccd_onl_hc <- make_results(ONL_CCD_metaDat_hc, 'CCD', 'ONL') %>% mutate(coh = 'x2')
 ccd_onl_lc <- make_results(ONL_CCD_metaDat_lc, 'CCD', 'ONL') %>% mutate(coh = 'x0.5')

 nt_lab_hc  <- make_results(LAB_NT_metaDat_hc, 'NT', 'LAB')   %>% mutate(coh = 'x2')
 nt_lab_lc  <- make_results(LAB_NT_metaDat_lc, 'NT', 'LAB')   %>% mutate(coh = 'x0.5')
 nt_onl_hc  <- make_results(ONL_NT_metaDat_hc, 'NT', 'ONL')   %>% mutate(coh = 'x2')
 nt_onl_lc  <- make_results(ONL_NT_metaDat_lc, 'NT', 'ONL')   %>% mutate(coh = 'x0.5')

meta_d_df_coh <- rbind(ccd_lab_hc, ccd_lab_lc, ccd_onl_hc, ccd_onl_lc,
                   nt_lab_hc, nt_lab_lc, nt_onl_hc, nt_onl_lc)
}

### Test of differences ----

m_rat_coh_test <- meta_d_df_coh %>%
  filter(str_detect(name, "Mratio"), !name %in% c('mu_logMratio', 'sigma_logMratio')) %>%
  mutate(ID = c(rep(1:11, 2), rep(1:10, 2), rep(1:7, 2), rep(1:84, 2)))

summary(lmerTest::lmer(mean ~ coh * group + (1|ID), data = m_rat_coh_test))

m_d1_coh_test <- meta_d_df_coh %>%
  filter(str_detect(name, "meta_d"), !name %in% c('mu_meta_d')) %>%
  mutate(ID = c(rep(1:11, 2), rep(1:10, 2), rep(1:7, 2), rep(1:84, 2)))

summary(lmerTest::lmer(mean ~ group + (1|ID),
           m_d1_coh_test %>%
             filter(coh == 'x0.5',
                    type == 'LAB')))

summary(lmerTest::lmer(mean ~ group * coh + (1|ID),
           m_d1_coh_test %>%
             filter(type == 'LAB')))

summary(lmerTest::lmer(mean ~ group + (1|ID),
           m_d1_coh_test %>%
             filter(coh == 'x0.5',
                    type == 'ONL')))

### Confidence spread ----

conf_spread <- data.frame(
  conf = c(6:1, -1:-6),
  cond = c(rep('lab', 72), rep('online', 72)),
  type = rep(c(rep('x2', 24), rep('Both', 24), rep('x0.5', 24)), 2),
  s1 = c(rowSums(LAB_CCD_metaDat_hc$S1mat), rowSums(LAB_NT_metaDat_hc$S1mat),
         rowSums(LAB_CCD_metaDat$S1mat),    rowSums(LAB_NT_metaDat$S1mat),
         rowSums(LAB_CCD_metaDat_lc$S1mat), rowSums(LAB_NT_metaDat_lc$S1mat),
         rowSums(ONL_CCD_metaDat_hc$S1mat), rowSums(ONL_NT_metaDat_hc$S1mat),
         rowSums(ONL_CCD_metaDat$S1mat),    rowSums(ONL_NT_metaDat$S1mat),
         rowSums(ONL_CCD_metaDat_lc$S1mat), rowSums(ONL_NT_metaDat_lc$S1mat)),
  s2 = c(rev(rowSums(LAB_CCD_metaDat_hc$S2mat)), rev(rowSums(LAB_NT_metaDat_hc$S2mat)),
         rev(rowSums(LAB_CCD_metaDat$S2mat)), rev(rowSums(LAB_NT_metaDat$S2mat)),
         rev(rowSums(LAB_NT_metaDat_lc$S2mat)), rev(rowSums(LAB_NT_metaDat_lc$S2mat)),
         rev(rowSums(ONL_CCD_metaDat_hc$S2mat)), rev(rowSums(ONL_NT_metaDat_hc$S2mat)),
         rev(rowSums(ONL_CCD_metaDat$S2mat)), rev(rowSums(ONL_NT_metaDat$S2mat)),
         rev(rowSums(ONL_NT_metaDat_lc$S2mat)), rev(rowSums(ONL_NT_metaDat_lc$S2mat)))
) %>%
  pivot_longer(s1:s2, names_to = 'Dec', values_to = 'Count') %>%
  mutate(group = c(rep('CCD', 144), rep('NT', 144))) %>%
  group_by(type, Dec, cond, group) %>%
  mutate(type = factor(type, levels = c('x0.5', 'Both', 'x2')),
         Count = Count/sum(Count)) %>%
  ggplot(aes(conf, Count, fill = group))+
    geom_col(position = 'dodge')+
    labs(x = 'Confidence Bin (1=low, 6=high)', y = 'Normalised Count')+
    scale_x_continuous(breaks = c(6, 1, -1, -6))+
    facet_wrap(~type)+
    theme_bw(base_size = 24)+
    theme(legend.position = c(0.10, 0.8),
          legend.background = element_blank(),
          legend.title = element_blank(),
          panel.grid = element_blank())

conf_spread

### contrast effects ------

contrast <- rbind(meta_d_df %>%
                    mutate(coh='both'),
                  meta_d_df_coh) %>%
  mutate(coh = factor(coh, levels = c('x0.5', 'both', 'x2')))

contrast_plot <- ggplot(contrast %>%
         filter(str_detect(name, "meta_d"),
                name!='mu_meta_d',
                type %in% c('ONL', 'LAB')),
       aes(x = coh, y = abs(mean), fill = group))+
  geom_boxplot(width = 0.4, outliers = F)+
  #stat_compare_means(hide.ns = T, label = 'p.signif')+
  facet_wrap(~type, scales = 'free_y')+
  labs(x = 'Coherence', y = "meta-d'")+

ggplot(contrast %>%
         filter(str_detect(name, "Mratio"),
                !name%in%c('mu_logMratio', 'sigma_logMratio'),
                type %in% c('ONL', 'LAB')),
       aes(x = coh, y = mean, fill = group))+
  coord_cartesian(ylim = c(0, 3))+
  geom_boxplot(width = 0.4, outliers = F)+
  #stat_compare_means(hide.ns = T, label = 'p.signif', label.y = 2.8)+
  labs(x = 'Coherence', y = "m-ratio")+
  facet_wrap(~type)&
  theme_bw(base_size = 24)&
  theme(legend.position = 'none',
        panel.grid = element_blank())

conf_spread/contrast_plot

# H-meta D permutation ---------------------------------------------------------

# Extract the data for the permutation analysis
mu_log_onl_NT             <- nt_onl  %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_onl_CCD            <- ccd_onl %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_onl_NT_d1          <- nt_onl  %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_onl_CCD_d1         <- ccd_onl %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_onl_NT_meta_d      <- nt_onl  %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))
mu_log_onl_CCD_meta_d     <- ccd_onl %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))

mu_log_lab_NT          <- nt_lab  %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_lab_CCD         <- ccd_lab %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_lab_NT_d1       <- nt_lab  %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_lab_CCD_d1      <- ccd_lab %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_lab_NT_meta_d   <- nt_lab  %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))
mu_log_lab_CCD_meta_d  <- ccd_lab %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))

mu_log_vr_NT              <- nt_b  %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_CCD             <- ccd_b %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_NT_d1           <- nt_b  %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_d1          <- ccd_b %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_NT_meta_d       <- nt_b  %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_meta_d      <- ccd_b %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))

mu_log_vr_NT_L            <- nt_l %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_CCD_L           <- ccd_l  %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_NT_d1_L         <- nt_l %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_d1_L        <- ccd_l  %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_NT_meta_d_L     <- nt_l %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_meta_d_L    <- ccd_l  %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))

mu_log_vr_NT_M            <- nt_m  %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_CCD_M           <- ccd_m   %>% filter(name == 'mu_logMratio') %>% mutate(mean = exp(mean))
mu_log_vr_NT_d1_M         <- nt_m  %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_d1_M        <- ccd_m   %>% filter(name == 'mu_d1') %>% mutate(mean = abs(mean))
mu_log_vr_NT_meta_d_M     <- nt_m  %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))
mu_log_vr_CCD_meta_d_M    <- ccd_m   %>% filter(name == 'mu_meta_d') %>% mutate(mean = abs(mean))

diffLabCCD         = as.numeric(mu_log_lab_NT[2] - mu_log_lab_CCD[2])
diffLabCCD_d1      = as.numeric(mu_log_lab_NT_d1[2] - mu_log_lab_CCD_d1[2])
diffLabCCD_meta_d  = as.numeric(mu_log_lab_NT_meta_d[2] - mu_log_lab_CCD_meta_d[2])

diffOnlCCD         = as.numeric(mu_log_onl_NT[2] - mu_log_onl_CCD[2])
diffOnlCCD_d1      = as.numeric(mu_log_onl_NT_d1[2] - mu_log_onl_CCD_d1[2])
diffOnlCCD_meta_d  = as.numeric(mu_log_onl_NT_meta_d[2] - mu_log_onl_CCD_meta_d[2])

diffVRCCD          = as.numeric(mu_log_vr_NT[2]  - mu_log_vr_CCD[2])
diffVRCCD_d1       = as.numeric(mu_log_vr_NT_d1[2]  - mu_log_vr_CCD_d1[2])
diffVRCCD_meta_d   = as.numeric(mu_log_vr_NT_meta_d[2]  - mu_log_vr_CCD_meta_d[2])

diffVRCCD_L        = as.numeric(mu_log_vr_NT_L[2]  - mu_log_vr_CCD_L[2])
diffVRCCD_d1_L     = as.numeric(mu_log_vr_NT_d1_L[2]  - mu_log_vr_CCD_d1_L[2])
diffVRCCD_meta_d_L = as.numeric(mu_log_vr_NT_meta_d_L[2]  - mu_log_vr_CCD_meta_d_L[2])

diffVRCCD_M        = as.numeric(mu_log_vr_NT_M[2]  - mu_log_vr_CCD_M[2])
diffVRCCD_d1_M     = as.numeric(mu_log_vr_NT_d1_M[2]  - mu_log_vr_CCD_d1_M[2])
diffVRCCD_meta_d_M = as.numeric(mu_log_vr_NT_meta_d_M[2]  - mu_log_vr_CCD_meta_d_M[2])

cat("The difference is mu", diffLabCCD, ' for Lab Vs CCD \n')
cat("The difference is d1", diffLabCCD_d1, ' for Lab Vs CCD \n')
cat("The difference is meta_d", diffLabCCD_meta_d, ' for Lab Vs CCD \n')

cat("The difference is ", diffOnlCCD, ' for Onl Vs CCD \n')

cat("The difference is mu ", diffVRCCD, ' for VR Vs CCD \n')
cat("The difference is d1", diffVRCCD_d1, ' for VR Vs CCD \n')
cat("The difference is meta_d", diffVRCCD_meta_d, ' for VR Vs CCD \n')

cat("The difference is mu", diffVRCCD_L, ' for VR Vs CCD \n')
cat("The difference is d1", diffVRCCD_d1_L, ' for VR Vs CCD \n')
cat("The difference is meta_d", diffVRCCD_meta_d_L, ' for VR Vs CCD \n')

cat("The difference is mu", diffVRCCD_M, ' for VR Vs CCD \n')
cat("The difference is d1", diffVRCCD_d1_M, ' for VR Vs CCD \n')
cat("The difference is meta_d", diffVRCCD_meta_d_M, ' for VR Vs CCD \n')

##Data for permutation -----

data_for_permute_lab_CCD     <- LAB_CCDRDKcheck %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)
data_for_permute_lab_NT      <- LAB_NTRDKcheck %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)

data_for_permute_onl_CCD     <- ONL_CCDRDKcheck %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)
data_for_permute_onl_NT      <- ONL_NTRDKcheck %>% dplyr::select(ID, Trial, conf, dotDirection, type, kmed, referenceSelection, group)

data_for_permute_vr_CCD      <- data_vr %>% filter(BinocularTrial == 1, Group == 'CCD')
data_for_permute_vr_NT       <- data_vr %>% filter(BinocularTrial == 1, Group == 'NT')

data_for_permute_vr_CCD_L    <- data_vr %>% filter(LateralizedTrial == 1, Group == 'CCD')
data_for_permute_vr_NT_L     <- data_vr %>% filter(LateralizedTrial == 1, Group == 'NT')

data_for_permute_vr_CCD_M    <- data_vr %>% filter(MonocularTrial == 1, Group == 'CCD')
data_for_permute_vr_NT_M     <- data_vr %>% filter(MonocularTrial == 1, Group == 'NT')

## Run the permutation -----

# Run the permutation function; NB: this will take a while
setwd('HMeta-d-master/R/')
source('Function_metad_group.R')

nreps <- 500

#lab_CCDvsNTp <- permute_metad_group(data_for_permute_lab_CCD, data_for_permute_lab_NT, cores = 10, nreps = nreps);
#onl_CCDvsNTp <- permute_metad_group(data_for_permute_onl_CCD, data_for_permute_onl_NT, cores = 10, nreps = nreps);
#vr_CCDvsNTp  <- permute_metad_group_vr(data_for_permute_vr_CCD, data_for_permute_vr_NT,cores = 10, nreps = nreps);
#vr_CCDvsNTp_L<- permute_metad_group_vr(data_for_permute_vr_CCD_L, data_for_permute_vr_NT_L,cores = 10, nreps = nreps);
#vr_CCDvsNTp_M<- permute_metad_group_vr(data_for_permute_vr_CCD_M, data_for_permute_vr_NT_M,cores = 10, nreps = nreps);

#saveRDS(lab_CCDvsNTp,  'HMetaDFit/CCDvsNT_LAB.rdata')
#saveRDS(onl_CCDvsNTp,  'HMetaDFit/CCDvsNT_ONL.rdata')
#saveRDS(vr_CCDvsNTp,  'HMetaDFit/CCDvsNT_VR.rdata')
#saveRDS(vr_CCDvsNTp_L,'HMetaDFit/CCDvsNT_VR_L.rdata')
#saveRDS(vr_CCDvsNTp_M,'HMetaDFit/CCDvsNT_VR_M.rdata')

## Load the permutation -----

lab_CCDvsNTp     <- readRDS('HMetaDFit/CCDvsNT_LAB.rdata')
onl_CCDvsNTp     <- readRDS('HMetaDFit/CCDvsNT_ONL.rdata')
vr_CCDvsNTp      <- readRDS('HMetaDFit/CCDvsNT_VR.rdata')
vr_CCDvsNTp_L    <- readRDS('HMetaDFit/CCDvsNT_VR_L.rdata')
vr_CCDvsNTp_M    <- readRDS('HMetaDFit/CCDvsNT_VR_M.rdata')

## Extract differences -----

Labdiff.obt       <- round(diffLabCCD, digits = 2)
Labdiff_d1.obt    <- round(diffLabCCD_d1, digits = 2)
Labdiff_metad.obt <- round(diffLabCCD_meta_d, digits = 2)

#effect size and superiority:
perm_effect_stats(Labdiff.obt, lab_CCDvsNTp[,4])
perm_effect_stats(Labdiff_d1.obt, lab_CCDvsNTp[,7])
perm_effect_stats(Labdiff_metad.obt, lab_CCDvsNTp[,10])

Onldiff.obt       <- round(diffOnlCCD, digits = 2)
Onldiff_d1.obt    <- round(diffOnlCCD_d1, digits = 2)
Onldiff_metad.obt <- round(diffOnlCCD_meta_d, digits = 2)

#effect size and superiority:
perm_effect_stats(Onldiff.obt, onl_CCDvsNTp[,4])
perm_effect_stats(Onldiff_d1.obt, onl_CCDvsNTp[,7])
perm_effect_stats(Onldiff_metad.obt, onl_CCDvsNTp[,10])

VRdiff.obt        <- round(diffVRCCD, digits = 2)
VRdiff_d1.obt     <- round(diffVRCCD_d1, digits = 2)
VRdiff_metad.obt  <- round(diffVRCCD_meta_d, digits = 2)

perm_effect_stats(VRdiff.obt, vr_CCDvsNTp[,4])
perm_effect_stats(VRdiff_d1.obt, vr_CCDvsNTp[,7])
perm_effect_stats(VRdiff_metad.obt, vr_CCDvsNTp[,10])

VRdiff.obt_L      <- round(diffVRCCD_L, digits = 2)
VRdiff_d1.obt_L   <- round(diffVRCCD_d1_L, digits = 2)
VRdiff_metad.obt_L<- round(diffVRCCD_meta_d_L, digits = 2)

perm_effect_stats(VRdiff.obt_L, vr_CCDvsNTp_L[,4])
perm_effect_stats(VRdiff_d1.obt_L, vr_CCDvsNTp_L[,7])
perm_effect_stats(VRdiff_metad.obt_L, vr_CCDvsNTp_L[,10])

VRdiff.obt_M      <- round(diffVRCCD_M, digits = 2)
VRdiff_d1.obt_M   <- round(diffVRCCD_d1_M, digits = 2)
VRdiff_metad.obt_M<- round(diffVRCCD_meta_d_M, digits = 2)

perm_effect_stats(VRdiff.obt_M, vr_CCDvsNTp_M[,4])
perm_effect_stats(VRdiff_d1.obt_M, vr_CCDvsNTp_M[,7])
perm_effect_stats(VRdiff_metad.obt_M, vr_CCDvsNTp_M[,10])

## Plot outcomes -----

p_p1 <- lab_CCDvsNTp %>%
  as.data.frame() %>%
  rename(Diff = 4) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = diffLabCCD + 0.07, y = 10, label = Labdiff.obt, fill = '#475B5A', colour ='white') +
  geom_vline(xintercept = c(diffLabCCD), color = c('#475B5A'), size = 1.5) +
  geom_text(x = -0.4, y = 15.0, label = paste('p = ', probLab), check_overlap = T, fontface = 'bold') +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        plot.margin = margin(1,1,1,1, unit = 'cm'))

p_p1

p_p2 <- onl_CCDvsNTp %>%
  as.data.frame() %>%
  rename(Diff = 4) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = diffOnlCCD + 0.07, y = 10, label = Onldiff.obt, fill = '#81D6E3', colour = 'white') +
  geom_vline(xintercept = c(diffOnlCCD), color = c('#81D6E3'), size = 1.5) +
  geom_text(x = -0.6, y = 10.0, label = paste('p = ', probOnl), check_overlap = T, fontface = 'bold') +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = 'cm'))

p_p2

(p_p2|p_p1)   # Makes level tags all uppercase

p_p_b <- vr_CCDvsNTp %>%
  as.data.frame() %>%
  rename(Diff = 4) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = diffVRCCD + 0.15, y = 7.5, label = VRdiff.obt, fill = colpal_vr[2], colour = 'white')+
  geom_vline(xintercept = c(diffVRCCD), color = c(colpal_vr[2]), size = 1.5) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = 'cm'))

p_p_b

p_p_l <- vr_CCDvsNTp_L %>%
  as.data.frame() %>%
  rename(Diff = 4) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = diffVRCCD_L + 0.15, y = 7.5, label = VRdiff.obt_L, fill = colpal_vr[2], colour = 'white')+
  geom_vline(xintercept = c(diffVRCCD_L), color = c(colpal_vr[2]), size = 1.5) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = 'cm'))

p_p_l

p_p_m <- vr_CCDvsNTp_M %>%
  as.data.frame() %>%
  rename(Diff = 4) %>%
  ggplot(aes(Diff)) +
  geom_histogram(binwidth = 0.01, position = 'dodge') +
  labs(
       x = expression(paste("Randomised ",mu['ratio'])),
       y = 'Sample Count') +
  geom_label(x = diffVRCCD_M + 0.15, y = 7.5, label = VRdiff.obt_M, fill = colpal_vr[2], colour = 'white')+
  geom_vline(xintercept = c(diffVRCCD_M), color = c(colpal_vr[2]), size = 1.5) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = 'cm'))

p_p_m


## Summary Heatmap of outcomes ---------------------------------------------

library(dplyr)
library(ggplot2)

# Example: build the table manually from your stored values
heat_data <- tribble(
  ~experiment, ~measure, ~diff, ~U3,
  "Lab",   "m-ratio",   Labdiff.obt,       perm_effect_stats(Labdiff.obt, lab_CCDvsNTp[,4])$U,
  "Lab",   "d'",    Labdiff_d1.obt,    perm_effect_stats(Labdiff_d1.obt, lab_CCDvsNTp[,7])$U,
  "Lab",   "meta-d'",Labdiff_metad.obt, perm_effect_stats(Labdiff_metad.obt, lab_CCDvsNTp[,10])$U,

  "Online", "m-ratio",   Onldiff.obt,       perm_effect_stats(Onldiff.obt, onl_CCDvsNTp[,4])$U,
  "Online", "d'",    Onldiff_d1.obt,    perm_effect_stats(Onldiff_d1.obt, onl_CCDvsNTp[,7])$U,
  "Online", "meta-d'",Onldiff_metad.obt, perm_effect_stats(Onldiff_metad.obt, onl_CCDvsNTp[,10])$U,

  "VR (BIN)",     "m-ratio",   VRdiff.obt,        perm_effect_stats(VRdiff.obt, vr_CCDvsNTp[,4])$U,
  "VR (BIN)",     "d'",    VRdiff_d1.obt,     perm_effect_stats(VRdiff_d1.obt, vr_CCDvsNTp[,7])$U,
  "VR (BIN)",     "meta-d'",VRdiff_metad.obt,  perm_effect_stats(VRdiff_metad.obt, vr_CCDvsNTp[,10])$U,

  "VR (LAT)",   "m-ratio",   VRdiff.obt_L,      perm_effect_stats(VRdiff.obt_L, vr_CCDvsNTp_L[,4])$U,
  "VR (LAT)",   "d'",    VRdiff_d1.obt_L,   perm_effect_stats(VRdiff_d1.obt_L, vr_CCDvsNTp_L[,7])$U,
  "VR (LAT)",   "meta-d'",VRdiff_metad.obt_L,perm_effect_stats(VRdiff_metad.obt_L, vr_CCDvsNTp_L[,10])$U,

  "VR (MON)",   "m-ratio",   VRdiff.obt_M,      perm_effect_stats(VRdiff.obt_M, vr_CCDvsNTp_M[,4])$U,
  "VR (MON)",   "d'",    VRdiff_d1.obt_M,   perm_effect_stats(VRdiff_d1.obt_M, vr_CCDvsNTp_M[,7])$U,
  "VR (MON)",   "meta-d'",VRdiff_metad.obt_M,perm_effect_stats(VRdiff_metad.obt_M, vr_CCDvsNTp_M[,10])$U
)

# Apply the white-out rule
heat_data <- heat_data %>%
  group_by(measure) %>%
  mutate(
    scaled_diff = (diff - min(diff, na.rm = TRUE)) /
                  (max(diff, na.rm = TRUE) - min(diff, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(fill_value = ifelse(U3 < 0.90, NA, scaled_diff),
         measure = factor(measure, levels = c('m-ratio', "d'", "meta-d'")),
         show_tile  = U3 >= 0.90,
         fill_value = ifelse(show_tile, scaled_diff, NA),
         label_txt  = ifelse(show_tile, sprintf("%.2f", diff), ""),  # print original diff
         text_col   = ifelse(!show_tile, "black",
                        ifelse(fill_value <= 0.25 | fill_value >= 0.75, "white", "black"))
         )

ggplot(heat_data, aes(x = measure, y = fct_rev(experiment), fill = fill_value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label_txt, color = text_col), fontface = "bold", size = 4) +
  scale_fill_gradient(
    limits = c(0, 1), low = 'lightblue', 'darkblue',
    na.value = "white",
  ) +
  scale_color_identity() +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')


## Paired stats ------------------------------------------------------------
nreps <- 500
# Within CCD
vr_CCDp_B_L  <- permute_metad_group_vr_within(data_for_permute_vr_CCD, data_for_permute_vr_CCD_L, cores = 10, nreps = nreps);
vr_CCDp_B_M  <- permute_metad_group_vr_within(data_for_permute_vr_CCD, data_for_permute_vr_CCD_M, cores = 10, nreps = nreps);
vr_CCDp_M_L  <- permute_metad_group_vr_within(data_for_permute_vr_CCD_L, data_for_permute_vr_CCD_M, cores = 10, nreps = nreps);

diffLabCCD_B_L     = as.numeric(mu_log_vr_CCD[2] - mu_log_vr_CCD_L[2])
diffLabCCD_B_L_d1  = as.numeric(mu_log_vr_CCD_d1[2] - mu_log_vr_CCD_d1_L[2])
diffLabCCD_B_L_md  = as.numeric(mu_log_vr_CCD_meta_d[2] - mu_log_vr_CCD_meta_d_L[2])

perm_effect_stats(length(vr_CCDp_B_L[,4][vr_CCDp_B_L[,4]   >= diffLabCCD_B_L])/nreps,    vr_CCDp_B_L[,4])
perm_effect_stats(length(vr_CCDp_B_L[,7][vr_CCDp_B_L[,7]   >= diffLabCCD_B_L_d1])/nreps, vr_CCDp_B_L[,7])
perm_effect_stats(length(vr_CCDp_B_L[,10][vr_CCDp_B_L[,10] >= diffLabCCD_B_L_md])/nreps, vr_CCDp_B_L[,10])

diffLabCCD_B_M     = as.numeric(mu_log_vr_CCD[2] - mu_log_vr_CCD_M[2])
diffLabCCD_B_M_d1  = as.numeric(mu_log_vr_CCD_d1[2] - mu_log_vr_CCD_d1_M[2])
diffLabCCD_B_M_md  = as.numeric(mu_log_vr_CCD_meta_d[2] - mu_log_vr_CCD_meta_d_M[2])

perm_effect_stats(length(vr_CCDp_B_M[,4][vr_CCDp_B_M[,4]   >= diffLabCCD_B_M])/nreps,    vr_CCDp_B_M[,4])
perm_effect_stats(length(vr_CCDp_B_M[,7][vr_CCDp_B_M[,7]   >= diffLabCCD_B_M_d1])/nreps, vr_CCDp_B_M[,7])
perm_effect_stats(length(vr_CCDp_B_M[,10][vr_CCDp_B_M[,10] >= diffLabCCD_B_M_md])/nreps, vr_CCDp_B_M[,10])

diffLabCCD_M_L     = as.numeric(mu_log_vr_CCD_L[2] - mu_log_vr_CCD_M[2])
diffLabCCD_M_L_d1  = as.numeric(mu_log_vr_CCD_d1_L[2] - mu_log_vr_CCD_d1_M[2])
diffLabCCD_M_L_md  = as.numeric(mu_log_vr_CCD_meta_d_L[2] - mu_log_vr_CCD_meta_d_M[2])

perm_effect_stats(length(vr_CCDp_B_L[,4][vr_CCDp_B_L[,4]   >= diffLabCCD_B_L])/nreps,    vr_CCDp_B_L[,4])
perm_effect_stats(length(vr_CCDp_B_L[,7][vr_CCDp_B_L[,7]   >= diffLabCCD_B_L_d1])/nreps, vr_CCDp_B_L[,7])
perm_effect_stats(length(vr_CCDp_B_L[,10][vr_CCDp_B_L[,10] >= diffLabCCD_B_L_md])/nreps, vr_CCDp_B_L[,10])

# Within NT
vr_NTp_B_L  <- permute_metad_group_vr_within(data_for_permute_vr_NT, data_for_permute_vr_NT_L,cores = 10, nreps = nreps);
vr_NTp_B_M  <- permute_metad_group_vr_within(data_for_permute_vr_NT, data_for_permute_vr_NT_M,cores = 10, nreps = nreps);
vr_NTp_M_L  <- permute_metad_group_vr_within(data_for_permute_vr_NT_L, data_for_permute_vr_NT_M,cores = 10, nreps = nreps);

diffLabNT_B_L  = as.numeric(mu_log_vr_NT[2] - mu_log_vr_NT_L[2])
diffLabNT_B_M  = as.numeric(mu_log_vr_NT[2] - mu_log_vr_NT_M[2])
diffLabNT_L_M  = as.numeric(mu_log_vr_NT_L[2] - mu_log_vr_NT_M[2])

perm_effect_stats(length(vr_NTp_B_L[,4][vr_NTp_B_L[,4]   >= diffLabNT_B_L])/nreps, vr_NTp_B_L[,4])
perm_effect_stats(length(vr_NTp_B_M[,7][vr_NTp_B_M[,7]   >= diffLabNT_B_M])/nreps, vr_NTp_B_M[,7])
perm_effect_stats(length(vr_NTp_M_L[,10][vr_NTp_M_L[,10] >= diffLabNT_L_M])/nreps, vr_NTp_M_L[,10])

# Supplement ----

## Calibration Check ----

### corrs with confidence ----

check_coh_lm <- checkBothRDK %>%
  group_by(ID) %>%
  filter(type == 'main', conf!=0) %>%
  mutate(kmed_av = ifelse(kmed=='kmed x 0.5', coherence*2, coherence/2),
         conf = mean(conf)) %>%
  dplyr::select(kmed_av, conf, ID, group) %>%
  distinct()

ggplot(aes(conf, kmed_av, colour=group), data = check_coh_lm) +
  geom_point()+
  geom_smooth(method = 'lm') +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_color_manual(values = colpal)+
  stat_cor(show.legend = F, method = 'spearman')+
  stat_cor(show.legend = F, aes(conf, kmed_av), colour = 'black', data = check_coh_lm, label.y.npc = 0.50, method = 'spearman')+
  theme_bw(base_size=20)+
  theme(legend.position = 'left',
        legend.title = element_blank()) |

ggplot(checkBothRDK %>%
            group_by(ID) %>%
            filter(type == 'main', conf!=0, kmed != 'kmed', coherence <= 0.5),
       aes(factor(conf), coherence, fill=group)) +
  geom_col(position = 'dodge')+
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_fill_manual(values = colpal)+
  scale_colour_manual(values = colpal)+
  facet_wrap(~kmed)+
  stat_cor(show.legend = F,
           method = 'spearman',
           aes(conf, coherence, colour = group))+
  theme_bw(base_size=20)+
  theme(legend.position = 'none')

### plots ----

# Plot group mean coherence across the calibration trials
check_coh_RDK <- rbind(us_CCDRDKcheck, au_CCDRDKcheck, online_NTRDKcheck, lab_NTRDKcheck) %>%
  mutate(prestype = NA,
         vr = 0) %>%
  dplyr::select(ID, Trial, coherence, group, prestype, type, vr) %>%
  rbind(.,
           data_vr %>%
             rename(ID = PID,
                    Trial = Trial,
                    group = Group,
                    coherence = ActiveCoherence,
                    prestype = PresentationType,
                    type = TrialType,
                    ) %>%
          mutate(vr=1) %>%
          dplyr::select(ID, Trial, coherence, group, prestype, type, vr))

CCDcheck2<-ggExtra::ggMarginal(
  ggplot(check_coh_RDK %>%
           filter(type=='calibration', vr==0))+
    stat_summary(aes(Trial, coherence, fill = group), geom = 'ribbon', alpha = 0.1, colour = NA)+
    stat_summary(aes(Trial, coherence, colour = group), geom = 'line')+
    geom_point(data = check_coh_RDK %>%
                 filter(type=='calibration', Trial == 118),
               aes(x=120, y=coherence, group = group, colour = group),
               shape = 22, size = 0, show.legend = F)+
    scale_color_manual(values = colpal)+
    scale_fill_manual(values = colpal)+
    labs(x = 'Trial', y = 'Coherence')+
    coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 118))+
    scale_x_continuous(expand = c(0,0))+
    theme_minimal(base_size=20)+
    theme(legend.position = 'top',
          panel.grid.minor = element_blank(),
          legend.title = element_blank()),
  margins = 'y', groupFill = T)


bino_check <-   ggExtra::ggMarginal(
  ggplot(check_coh_RDK %>%
           filter(type=='Calibration',
                  prestype=='Binocular',
                  vr==1) %>%
           group_by(ID) %>%
           mutate(Trial = 1:n()))+
    stat_summary(aes(Trial, coherence, fill = group), geom = 'ribbon', alpha = 0.1, colour = NA)+
    stat_summary(aes(Trial, coherence, colour = group), geom = 'line')+
    geom_point(data = check_coh_RDK %>%
                 filter(type=='Calibration', Trial == 24),
               aes(x=120, y=coherence, group = group, colour = group),
               shape = 22, size = 0, show.legend = F)+
    scale_color_manual(values = colpal_vr, name = 'Bino.')+
    scale_fill_manual(values = colpal_vr , name = 'Bino.')+
    labs(x = 'Trial', y = 'Coherence')+
    coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 24))+
    scale_x_continuous(expand = c(0,0))+
    theme_minimal(base_size=20)+
    theme(legend.position = 'top',
          panel.grid.minor = element_blank()),
  margins = 'y', groupFill = T)

lat_check <-  ggExtra::ggMarginal(
  ggplot(check_coh_RDK %>%
           filter(type=='Calibration',
                  prestype=='Lateralized',
                  vr==1) %>%
           group_by(ID) %>%
           mutate(Trial = 1:n()))+
    stat_summary(aes(Trial, coherence, fill = group), geom = 'ribbon', alpha = 0.1, colour = NA)+
    stat_summary(aes(Trial, coherence, colour = group), geom = 'line')+
    geom_point(data = check_coh_RDK %>%
                 filter(type=='Calibration', Trial == 48),
               aes(x=120, y=coherence, group = group, colour = group),
               shape = 22, size = 0, show.legend = F)+
    scale_color_manual(values = colpal_vr, name = 'Lat.')+
    scale_fill_manual(values = colpal_vr , name = 'Lat.')+
    labs(x = 'Trial', y = 'Coherence')+
    coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 48))+
    scale_x_continuous(expand = c(0,0))+
    theme_minimal(base_size=20)+
    theme(legend.position = 'top',
          panel.grid.minor = element_blank()),
  margins = 'y', groupFill = T)

mono_check <-  ggExtra::ggMarginal(
  ggplot(check_coh_RDK %>%
           filter(type=='Calibration',
                  prestype=='Monocular',
                  vr==1) %>%
           group_by(ID) %>%
           mutate(Trial = 1:n()))+
    stat_summary(aes(Trial, coherence, fill = group), geom = 'ribbon', alpha = 0.1, colour = NA)+
    stat_summary(aes(Trial, coherence, colour = group), geom = 'line')+
    geom_point(data = check_coh_RDK %>%
                 filter(type=='Calibration', Trial == 48),
               aes(x=120, y=coherence, group = group, colour = group),
               shape = 22, size = 0, show.legend = F)+
    scale_color_manual(values = colpal_vr, name = 'Mono.')+
    scale_fill_manual(values = colpal_vr , name = 'Mono.')+
    labs(x = 'Trial', y = 'Coherence')+
    coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 48))+
    scale_x_continuous(expand = c(0,0))+
    theme_minimal(base_size=20)+
    theme(legend.position = 'top',
          panel.grid.minor = element_blank()),
  margins = 'y', groupFill = T)

coh_vr_check_full <- ggarrange(bino_check,lat_check,mono_check, nrow=1)

CCDcheck2
coh_vr_check_full

## Raw plot of Exp 1 & 2 group diffs ----

tidy_acc <- tidyplot(av_cor %>%
           mutate(
              kmed = ifelse(kmed=='kmed x 0.5', 'x0.5', 'x2.0'),
              kmed = as.character(kmed),
              group = ifelse(group=='CCD (MRI)', 'CCD\n(MRI)',
                             ifelse(group=='NT (Online)', 'NT\n(Online)',
                                    ifelse(group == 'CCD (Online)', 'CCD\n(Online)', 'NT\n(MRI)')))),
                  x = group, y = av_cor, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) %>%
  add_sd_errorbar() %>%
  add_mean_dot(size=2) %>%
  adjust_size(width = 75, height = 50) %>%
  adjust_x_axis_title(title = '') %>%
  adjust_y_axis_title(title = 'p(Correct)') %>%
  adjust_font(fontsize = 15) %>%
  adjust_legend_title(title = 'Coh.') %>%
  adjust_legend_position(position = 'top')+
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 1.04) +
  coord_cartesian(clip = "off")

tidy_conf <- tidyplot(av_co %>%
           mutate(
              kmed = ifelse(kmed=='kmed x 0.5', 'x0.5', 'x2.0'),
              kmed = as.character(kmed),
              group = ifelse(group=='CCD (MRI)', 'CCD\n(MRI)',
                             ifelse(group=='NT (Online)', 'NT\n(Online)',
                                    ifelse(group == 'CCD (Online)', 'CCD\n(Online)', 'NT\n(MRI)'))),
              ID = factor(ID)),
                  x = group, y = av_co, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) %>%
  add_sd_errorbar() %>%
  add_mean_dot(size=2) %>%
  adjust_size(width = 75, height = 50) %>%
  adjust_x_axis_title(title = '') %>%
  adjust_y_axis_title(title = 'Confidence') %>%
  adjust_font(fontsize = 15) %>%
  adjust_legend_title(title = 'Coh.') %>%
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 102) +
  coord_cartesian(clip = "off")

tidy_acc | tidy_conf

## Raw plot of Exp 3 group diffs ----

tidyplot_exp3_co <- av_co_vr %>%
           mutate(
              kmed = as.character(kmed),
              Pres = ifelse(PresentationType=='Monocular', 'Mono.',
                             ifelse(PresentationType=='Lateralized', 'Lat.', 'Bino.')),
              Pres = factor(Pres),
              ID = factor(ID))

tidyplot_exp3_cor <- av_cor_vr %>%
           mutate(
              kmed = as.character(kmed),
              Pres = ifelse(PresentationType=='Monocular', 'Mono.',
                             ifelse(PresentationType=='Lateralized', 'Lat.', 'Bino.')),
              Pres = factor(Pres),
              ID = factor(ID))

Conf <- tidyplot(tidyplot_exp3 %>% filter(Pres=='Bino.'),
                  x = group, y = av_co, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'Confidence') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 0.75) +
  coord_cartesian(clip = "off") |

tidyplot(tidyplot_exp3 %>% filter(Pres=='Lat.'),
                  x = group, y = av_co, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'Confidence') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 0.75) +
  coord_cartesian(clip = "off") |

tidyplot(tidyplot_exp3 %>% filter(Pres=='Mono.'),
                  x = group, y = av_co, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'Confidence') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 0.75) +
  coord_cartesian(clip = "off")

Cor <-tidyplot(tidyplot_exp3_cor %>% filter(Pres=='Bino.'),
                  x = group, y = av_cor, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  add_title('Bino.') |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'p(Correct)') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 1.04) +
  coord_cartesian(clip = "off") |

tidyplot(tidyplot_exp3_cor %>% filter(Pres=='Lat.'),
                  x = group, y = av_cor, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  add_title('Lat.') |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'p(Correct)') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 1.04) +
  coord_cartesian(clip = "off") |

tidyplot(tidyplot_exp3_cor %>% filter(Pres=='Mono.'),
                  x = group, y = av_cor, colour = kmed) %>%
  add_data_points_beeswarm(alpha = 0.1, cex = 1) |>
  add_sd_errorbar() |>
  add_mean_dot(size=2) |>
  add_title('Mono.') |>
  adjust_size(width = 75, height = 50) |>
  adjust_x_axis_title(title = '') |>
  adjust_y_axis_title(title = 'p(Correct)') |>
  adjust_font(fontsize = 15) |>
  adjust_legend_title(title = 'Coh.') |>
  adjust_legend_position(position = 'top') +
  stat_compare_means(paired = T,
                     hide.ns = T,
                     size = 5,
                     show.legend = F,
                     label = 'p.signif',
                     label.y = 1.04) +
  coord_cartesian(clip = "off")

Cor/Conf
