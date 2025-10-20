# 13.10.25
# Ethan Ross

# In this script I will plot the metaprobe gauze DNA concentrations pre and post PCR

library(ggplot2)
library(dplyr)
library(stringr)
library(writexl)
library(readxl)

met_conc <- read_xlsx("metaprobe_DNA_and_PCR_record.xlsx",
                      sheet = 2)

str(met_conc)

met_conc_long <- met_conc %>%
  pivot_longer(
    cols = c(Gauze, Pooled_Tele02, Pooled_MiFishU),
    names_to = "stage",
    values_to = "concentration"
  )

met_conc_plot <- ggplot(met_conc_long, aes(x = stage, y = concentration, color = factor(Deployment))) +
  geom_point(shape = 19, size = 4, alpha = 0.5) +
  labs(
    x = "",
    y = "Concentration (ng/ul)",
    fill = "Deployment"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C"),
    labels = c(
      "1" = "1 (Day 1)",
      "2" = "2 (Day 1)",
      "3" = "3 (Day 2)",
      "4" = "4 (Day 2)"
    )) +
  scale_y_sqrt()

met_conc_plot


################################################################################

met_conc_plot <- ggplot(met_conc_long, aes(x = stage, y = concentration)) +
  geom_boxplot(
    fill = "#F0F0F0", 
    outlier.shape = NA,
    alpha = 0.7,
    color = "black"
  ) +
  geom_jitter(
    aes(color = factor(Deployment)), 
    shape = 19, 
    size = 4,
    alpha = 0.8, 
    width = 0.2
  ) +
  scale_color_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C"),
    labels = c(
      "1" = "1 (Day 1)",
      "2" = "2 (Day 1)",
      "3" = "3 (Day 2)",
      "4" = "4 (Day 2)"
    )) +
  
  scale_y_log10() +
  labs(
    x = "",
    y = "DNA concentration (ng/µl)",
    color = "Deployment"
  ) +
  scale_x_discrete(
    labels = c(
      "Gauze" = "Gauze",
      "Pooled_Tele02" = "Pooled\nAmplicon\nTele02",
      "Pooled_MiFishU" = "Pooled\nAmplicon\nMiFishU"
    )
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17.5),
    axis.text.y = element_text(size = 17.5),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "right"
  )

print(met_conc_plot)

ggsave(filename = "supplimentary/metaprobe_DNA_concentration.png", plot = met_conc_plot,
       device = "png", dpi = 600, units = "mm", height = 150, width = 200)
