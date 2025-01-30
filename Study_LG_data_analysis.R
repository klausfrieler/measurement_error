### Study 2: Anaysis empirical data from LongGold
library(tidyverse)

LG_dat <- read.csv("~/Documents/Conferences/2025_ICMPC/select_vars.csv",colClasses=c("p_id"="character"))
#LG_dat <- read.csv("c:/Users/klaus.frieler/Downloads/select_vars.csv",colClasses=c("p_id"="character"))

### Pre-process data:
### use only cases where BAT, MIQ, GMS.MT are not missing
LG_dat <- LG_dat[!is.na(LG_dat$BAT.score) & !is.na(LG_dat$MIQ.score) & !is.na(LG_dat$GMS.musical_training), ]

### remove all cases without a valid ID (p_id not starting with "0")
LG_dat <- LG_dat %>% filter(substr(p_id, 1, 1) == "0")
  
### use each p_id only once, selecting teh observation from teh latest test wave or highest age
LG_dat <- LG_dat %>% 
  group_by(p_id) %>% 
  mutate(last = test_wave == max(test_wave)) %>% 
  ungroup() %>% 
  filter(last) 
### demographics

mean(LG_dat$age)
sd(LG_dat$age)
table(LG_dat$gender) / sum(table(LG_dat$gender))

### Data analysis with and without correction; only MI correction for now
require(TSI)
source("c:/Users/klaus.frieler/MPIAE/projects/measurement_error/my_TSI.R")
m_naive <- lm(BAT.score ~ MIQ.score + GMS.musical_training,data=LG_dat)

mice.dfs = my_TSI(
  LG_dat %>% select(BAT.score, BAT.error, MIQ.score, MIQ.error, GMS.musical_training),
  os_names = c("BAT.score", "MIQ.score"),
  se_names = c("BAT.error", "MIQ.error"),
  metrics = "z",
  score_types = "ML",
  separated = T,
  mice_args = c(
    m = 5,
    maxit = 10,
    printFlag = F
  )
)

coefs <- pool(with(mice.dfs, lm(true_BAT.score ~ true_MIQ.score + GMS.musical_training))) %>%
  summary() %>%
  select(term, value = estimate, se = std.error) %>%
  mutate(term = c("b0", "b1", "b2"))


