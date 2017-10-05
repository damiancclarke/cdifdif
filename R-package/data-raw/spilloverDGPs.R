library(tidyverse)


spilloverDGP <- read_csv("data-raw/spilloverDGPs.csv", col_names = FALSE)
spilloverDGP <- set_names(spilloverDGP, c("id", "dist", "x", "time", "treat", "y1", "y2", "y3"))
spilloverDGP

spilloverDGP2 <- spilloverDGP %>%
  filter(treat == 1) %>%
  distinct(id) %>%
  mutate(trettreat = 1) %>%
  left_join(spilloverDGP, . ) %>%
  mutate(trettreat = ifelse(is.na(trettreat), "0", trettreat))

ggplot(spilloverDGP2) +
  geom_line(aes(time, y1, group = id, color = trettreat))


devtools::use_data(spilloverDGP, overwrite = TRUE)
