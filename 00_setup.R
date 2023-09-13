library(tidyverse)
library(here)
library(readxl)
library(patchwork)
library(tools)
library(glue)

my_theme <- theme_grey() +
  theme(axis.text = element_text(size = 18)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(strip.text = element_text(face = "bold", size = 15)) +
  theme(legend.position = "right")
theme_set(my_theme)
