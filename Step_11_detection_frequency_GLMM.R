# 13.10.25
# Ethan Ross

# I will use a GLMM to test if detection by BRUV results in higher detection frequency


library(lme4)
library(readxl)
library(dplyr)

analysis_data <- read_xlsx("detection_frequency_data/metaprobe_detection_frequency_data.xlsx")

analysis_data <- analysis_data %>%
  mutate(
    BRUV = factor(BRUV),
    day = factor(day),
    deployment = factor(deployment)
  )

str(analysis_data)

model <- glmer(
  detection ~ BRUV + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

summary(model)
plot(model)

null_model <- glmer(
  detection ~ + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

anova(model, null_model, test = "Chisq")

################################################################################

str(analysis_data)

model_2 <- glmer(
  detection ~ BRUV + (1 | day) + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

summary(model_2)
plot(model_2)

null_model_2 <- glmer(
  detection ~ + (1 | day) + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

anova(model_2, null_model_2, test = "Chisq")

################################################################################

# Test all models

anova(model_2, model, test = "Chisq")

# Model 1 has a lower AIC so I will stick with it
# I.e. accounting for Day does not add much explanatory data to the model

# FINAL


summary(model)


