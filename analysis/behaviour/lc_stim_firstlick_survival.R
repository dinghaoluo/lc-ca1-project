# Created on 30 July 2026
# an R port of the lick-time survival analysis done with lifelines in Python 
# @author: Dinghao Luo 


## imports
library(tidyverse)
library(broom)
library(survival)
library(survRM2)
library(cmprsk)


## params 
repo_root     <- 'Z:/Dinghao/code/lc-ca1-project'
# cohort        <- 'soma'  # we are going to focus on LC-soma stims; terminal gave weak results 
entry_s       <- 1
horizon_run_s <- 5
n_bootstrap   <- 5000
seed          <- 42
set.seed(seed)


## load trial table which has been processed and saved by the Python script 
trial_table <- file.path(
  repo_root,
  'data',
  'behaviour',
  'lc_stim_firstlick_survival',
  'soma_firstlick_survival_trials.csv'
)

trials <- read_csv(
  trial_table
) |> 
  mutate(
    condition = factor(
      condition,
      levels = c('control', 'stim')  # good habit to use factors for control-treatment level pairs; probably not necessary
    )
  )


## we can directly use the time_from_landmark_s here which is simply duration_s minus entry_s 

## KM fitting 
km_fit <- survfit(
  Surv(time_from_landmark_s, event) ~ condition,
  data = trials,
  conf.type = 'log-log'  # this is the default for KaplanMeierFitter but survfit() uses 'log'
)

print(km_fit)  # this already gives stim median >> cont median 


## log-rank 
logrank <- survdiff(
  Surv(time_from_landmark_s, event) ~ condition,
  data = trials,
  rho = 0
)
print(logrank)  # gives p=2e-04


## Cox PH models
cox_basic <- coxph(
  Surv(time_from_landmark_s, event) ~ is_stim,
  data = trials,
  ties = 'efron'
)
summary(cox_basic)

cox_basic_tidy <- tidy(
  cox_basic,
  exponentiate = TRUE,  # convert estimate to hazard ratio 
  conf.int = TRUE
)
cox_basic_tidy

## add clustering to Cox 
cox_clustered <- coxph(
  Surv(time_from_landmark_s, event) ~ 
    is_stim + cluster(session),
  data = trials,
  ties = 'efron'
)
summary(cox_clustered)

cox_clustered_tidy <- tidy(
  cox_clustered,
  exponentiate = TRUE,
  conf.int = TRUE 
)
cox_clustered_tidy

ph_check <- cox.zph(cox_basic)  # this gives p=0.027; expected though, since stim would have a strong effect early on
plot(ph_check)


## RMST 
tau <- horizon_run_s - entry_s 
rmst_fit <- rmst2(
  time = trials$time_from_landmark_s,
  status = trials$event,
  arm = trials$is_stim,
  tau = tau 
)

print(rmst_fit)
str(rmst_fit, max.level=2)

rmst_control <- rmst_fit$RMST.arm0$result['RMST', 'Est.']
rmst_stim    <- rmst_fit$RMST.arm1$result['RMST', 'Est.']
rmst_diff    <- rmst_fit$unadjusted.result['RMST (arm=1)-(arm=0)', 'Est.']

print(c(
  control     = rmst_control,
  stimulation = rmst_stim,
  difference  = rmst_diff
))


## session bootstrap
sessions <- unique(trials$session)

bootstrap_rmst <- map_dbl(
  seq_len(n_bootstrap),
  function(iteration) {
    sampled_sessions <- sample(
      sessions,
      size = length(sessions),
      replace = TRUE
    )
    
    sampled_trials <- map_dfr(
      sampled_sessions,
      function(selected_session) {
        trials |>
          filter(session == selected_session)
      }
    )
    
    fit <- rmst2(
      time = sampled_trials$time_from_landmark_s,
      status = sampled_trials$event,
      arm = sampled_trials$is_stim,
      tau = tau
    )
    
    unname(
      fit$unadjusted.result[
        'RMST (arm=1)-(arm=0)',
        'Est.'
      ]
    )
  }
)

rmst_interval <- quantile(
  bootstrap_rmst,
  probs = c(0.025, 0.975),
  names = FALSE
)
print(rmst_interval)  # lower = 0.13 and upper = 0.47


## competing risk 
cif_fit <- cuminc(
  ftime = trials$time_from_landmark_s,
  fstatus = trials$event_type,
  group = trials$condition,
  cencode = 0  # 0 is censored 
)

print(cif_fit)
print(names(cif_fit))
print(cif_fit$Tests)

cif_at_time <- function(curve, time) {
  index <- max(which(curve$time <= time))
  curve$est[index]
}

cif_control <- cif_at_time(
  cif_fit[['control 1']],
  tau
)

cif_stim <- cif_at_time(
  cif_fit[['stim 1']],
  tau
)

cif_difference <- cif_stim - cif_control

print(c(
  control     = cif_control,
  stimulation = cif_stim,
  difference  = cif_difference
))


## plotting 
# paths 
figure_dir <- file.path(
  repo_root,
  'figures',
  'behaviour',
  'lc_stim_firstlick_survival'
)
dir.create(figure_dir, showWarnings = FALSE)

km_data <- tidy(km_fit) |>
  mutate(
    run_time_s = time + entry_s,
    condition = str_remove(strata, '^condition=')
  )

km_plot <- ggplot(
  km_data,
  aes(
    x = run_time_s,
    y = estimate,
    colour = condition,
    fill = condition
  )
) + 
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.25,
    colour = NA
  ) +
  geom_step(linewidth = 0.5) + 
  scale_colour_manual(
    values = c(control = 'grey', stim = 'navyblue')
  ) + 
  scale_fill_manual(
    values = c(control = 'grey', stim = 'navyblue')
  ) + 
  coord_cartesian(
    xlim = c(1,5),
    ylim = c(0,1)
  ) + 
  labs(
    x = 'Time from run onset (s)',
    y = 'P(not licked yet)',
    colour = NULL,
    fill = NULL
  )

ggsave(
  filename = file.path(
    figure_dir,
    'soma_firstlick_survival_r_first5s.png'
  ),
  plot = km_plot,
  width = 5.6,
  height = 3.4,
  dpi = 300
)


cif_names <- setdiff(
  names(cif_fit),
  'Tests'
)

cif_data <- map_dfr(
  cif_names,
  function(curve_name) {
    curve <- cif_fit[[curve_name]]
    
    tibble(
      time = curve$time,
      estimate = curve$est,
      variance = curve$var,
      curve = curve_name
    )
  }
) |>
  separate(
    curve,
    into = c('condition', 'event_type'),
    sep = ' '
  ) |>
  mutate(
    event_type = as.integer(event_type),
    run_time_s = time + entry_s,
    lower = pmax(0, estimate - 1.96 * sqrt(variance)),
    upper = pmin(1, estimate + 1.96 * sqrt(variance))
  ) |>
  filter(event_type == 1)

cif_plot <- ggplot(
  cif_data,
  aes(
    x = run_time_s,
    y = estimate,
    colour = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18,
    colour = NA
  ) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(
    values = c(control = 'grey', stim = 'navyblue')
  ) +
  scale_fill_manual(
    values = c(control = 'grey', stim = 'navyblue')
  ) +
  coord_cartesian(
    xlim = c(1, 5),
    ylim = c(0, 1)
  ) +
  labs(
    x = 'Time from run onset (s)',
    y = 'Cumulative first licks',
    colour = NULL,
    fill = NULL
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figure_dir,
    'soma_firstlick_cif_r_first5s.png'
  ),
  plot = cif_plot,
  width = 4.2,
  height = 3.4,
  dpi = 300
)