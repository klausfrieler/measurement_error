### Study 2: Anaysis empirical data from LongGold
library(tidyverse)

<<<<<<< HEAD
LG_dat <- read.csv("select_vars.csv",colClasses=c("p_id"="character"))
=======
LG_dat <- read.csv("~/Documents/Conferences/2025_ICMPC/select_vars.csv",colClasses=c("p_id"="character"))
#LG_dat <- read.csv("c:/Users/klaus.frieler/Downloads/select_vars.csv",colClasses=c("p_id"="character"))
>>>>>>> d3c50bb6625f485a2664085a784a07f03f5515ac

### Pre-process data:
### use only cases where BAT, MIQ, GMS.MT are not missing
LG_dat <- LG_dat[!is.na(LG_dat$BAT.score) & !is.na(LG_dat$MIQ.score) & !is.na(LG_dat$GMS.musical_training), ]

### remove all cases without a valid ID (p_id not starting with "0")
LG_dat <- LG_dat %>% filter(substr(p_id, 1, 1) == "0")
<<<<<<< HEAD

### use each p_id only once, selecting the observation from the latest test wave or highest age
=======
  
### use each p_id only once, selecting teh observation from teh latest test wave or highest age
>>>>>>> d3c50bb6625f485a2664085a784a07f03f5515ac
LG_dat <- LG_dat %>% 
  group_by(p_id) %>% 
  mutate(last = test_wave == max(test_wave)) %>% 
  ungroup() %>% 
  filter(last) 
<<<<<<< HEAD

=======
>>>>>>> d3c50bb6625f485a2664085a784a07f03f5515ac
### demographics

mean(LG_dat$age)
sd(LG_dat$age)
table(LG_dat$gender) / sum(table(LG_dat$gender))

### Data analysis with and without correction; only MI correction for now
require(TSI)
<<<<<<< HEAD
source("my_TSI.R")
require(tidyverse)

=======
source("c:/Users/klaus.frieler/MPIAE/projects/measurement_error/my_TSI.R")
>>>>>>> d3c50bb6625f485a2664085a784a07f03f5515ac
m_naive <- lm(BAT.score ~ MIQ.score + GMS.musical_training,data=LG_dat)
coefs_naive <- coef(m_naive)

<<<<<<< HEAD
#inv_err <- 1 / LG_dat$MIQ.error ^ 2
inv_err <- (1 / LG_dat$MIQ.error ^ 2) * (1 / LG_dat$BAT.error ^ 2) #* (1 / LG_dat$GMS.Musical_training.error ^ 2)
coefs_w <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = LG_dat)) %>% 
  select(term, value = estimate, se = std.error) %>% 
  mutate(term = c("b0", "b1", "b2"))
coefs_w

mice.dfs = my_TSI(
  LG_dat %>% select(-c(p_id, test_year,test_wave)),
=======
mice.dfs = my_TSI(
  LG_dat %>% select(BAT.score, BAT.error, MIQ.score, MIQ.error, GMS.musical_training),
>>>>>>> d3c50bb6625f485a2664085a784a07f03f5515ac
  os_names = c("BAT.score", "MIQ.score"),
  se_names = c("BAT.error", "MIQ.error"),
  metrics = "z",
  score_types = "ML",
  separated = T,
  mice_args = c(
    m = 20,
    maxit = 10,
    printFlag = F
  )
)

coefs_MI <- pool(with(mice.dfs, lm(true_BAT.score ~ true_MIQ.score + GMS.musical_training))) %>%
  summary() %>%
  select(term, value = estimate, se = std.error) %>%
  mutate(term = c("b0", "b1", "b2"))

coefs_MI
