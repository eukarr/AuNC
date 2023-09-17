source("00_setup.R")
source("01_functions.R")
# for subscript in ggplot
library(ggtext)
# ANOVA
library(car)


if (!dir.exists(here("plots"))) dir.create(here("plots"))

if (!dir.exists(here("plots", "mendcomm"))) dir.create(here("plots", "mendcomm"))
if (!dir.exists(here("plots", "analytica"))) dir.create(here("plots", "analytica"))



#### N3 series ####
samples_azide_ratio <- read_csv(here('data', 'N3', 'N3.csv'))


# import quantum yield fluorescence spectra

ratio_fluoro <- here('data', 'N3', '2023-07-14', 'fluoro') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*-", "", file))) %>%
  left_join(samples_azide_ratio, by = "file") %>%
  mutate(type = factor(type, levels = c("au", "ref"))) %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = glue("azide:gold {n3_au}"))



(fig_1a_mc <- ggplot(data = ratio_fluoro %>% filter(type != 'ref'), 
       aes(x = wavelength, y = intensity)) +
  geom_line(aes(color = factor(n3_au))) +
  lims(x = c(400, 800)) +
  labs(x = "Wavelength, nm", 
       y = "Intensity, a.u.", 
       color = "N<sub>3</sub> / Au") +
  theme(legend.title = element_markdown()))

ggsave(plot = fig_1a_mc,
       filename = here("plots", "mendcomm", "fig_1a.tiff"), 
       dpi = 300)



ratio_abs <- here('data', 'N3', '2023-07-14', 'abs') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*-", "", file))) %>%
  left_join(samples_azide_ratio, by = "file") %>%
  mutate(type = factor(type, levels = c("au", "ref"))) %>%
  filter(absorbance <= 1.9) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = glue("azide:gold {n3_au}"))


(fig_1c_mc <- ggplot(data = ratio_abs %>% filter(type != 'ref'), aes(x = wavelength, y = absorbance)) +
  geom_line(aes(color = factor(n3_au))) +
  lims(x = c(280, 500), y = c(0, 1.25)) +
  labs(x = "Wavelength, nm", 
       y = "Absorbance", 
       color = "N<sub>3</sub> / Au") +
    theme(legend.title = element_markdown()))

ggsave(plot = fig_1c_mc,
       filename = here("plots", "mendcomm", "fig_1c.tiff"), 
       dpi = 300)



ratio_abs365 <- ratio_abs %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(file) %>%
  summarize(absorbance = median(absorbance))

ratio_qy <- ratio_fluoro %>%
  left_join(ratio_abs365, by = "file")

ratio_integral <- ratio_qy %>%
  filter(wavelength > 400) %>%
  group_by(file, n3_au, type) %>%
  summarize(integral = sum(intensity / absorbance)) %>%
  ungroup()

ratio_ref_integral <- ratio_integral %>%
  filter(str_detect(type, "ref")) %>%
  select(type, integral) %>%
  rename(ref_integral = integral) %>%
  mutate(type = 'au')

ratio_qy <- ratio_integral %>%
  left_join(ratio_ref_integral, keep = FALSE, by = "type") %>%
  mutate(qy = integral / ref_integral * 0.54) %>%
  filter(!str_detect(type, "ref"))



(fig_1d_mc <- ggplot(data = ratio_qy) +
  geom_col(aes(x = factor(n3_au), y = qy * 100)) +
  labs(x = "N<sub>3</sub> / Au", y = "QY, %") +
  theme(axis.title.x = element_markdown(),
        axis.ticks.x = element_blank()))

ggsave(plot = fig_1d_mc,
       filename = here("plots", "mendcomm", "fig_1d.tiff"), 
       dpi = 300)


#### control experiments from factorial series ####
samples_factorial <- read_csv(here('data', 'factorial', 'factorial_AuNC.csv'))


# import absorbance spectra for quantum yield measurement
control_abs <- here('data', 'factorial', '2023-08-10', 'abs') %>%
  list.files(pattern = "^Absorbance.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type == 'control') %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = paste(ifelse(AMP_ratio, 'AMP', ' '),
                       ifelse(Citr_ratio, 'citr', ' '),
                       ifelse(N3_ratio, 'N<sub>3</sub>', ' '), sep = "+") %>%
           str_remove_all(fixed(" +")) %>%
           str_remove_all(fixed("+ ")))


(fig_2a_mc <- ggplot(data = control_abs, aes(x = wavelength, y = absorbance)) +
  geom_line(aes(color = label)) +
  lims(x = c(280, 500), y = c(0, 2.5)) +
  labs(x = "Wavelength, nm", 
       y = "Absorbance", 
       color = "Components") +
  theme(legend.text = element_markdown()))

ggsave(plot = fig_2a_mc,
       filename = here("plots", "mendcomm", "fig_2a.tiff"), 
       dpi = 300)


control_fluoro <- here('data', 'factorial', '2023-08-10', 'diluted_fluoro') %>%
  list.files(pattern = "^Subt2_.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type == 'control') %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = paste(ifelse(AMP_ratio, 'AMP', ' '),
                       ifelse(Citr_ratio, 'citr', ' '),
                       ifelse(N3_ratio, 'N<sub>3</sub>', ' '), sep = "+") %>%
           str_remove_all(fixed(" +")) %>%
           str_remove_all(fixed("+ ")))

(fig_2b_mc <- ggplot(data = control_fluoro, aes(x = wavelength, y = intensity)) +
    geom_line(aes(color = label)) +
    lims(x = c(400, 800)) +
    labs(x = "Wavelength, nm", 
         y = "Intensity, a.u.", 
         color = "Components") +
    theme(legend.text = element_markdown(),
          axis.text.y = element_blank()))

ggsave(plot = fig_2b_mc,
       filename = here("plots", "mendcomm", "fig_2b.tiff"), 
       dpi = 300)


control_fluoro_stored <- here('data', 'factorial', '2023-08-28', 'fluoro') %>%
  list.files(pattern = "^Subt2_.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type == 'control') %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = paste(ifelse(AMP_ratio, 'AMP', ' '),
                       ifelse(Citr_ratio, 'citr', ' '),
                       ifelse(N3_ratio, 'N<sub>3</sub>', ' '), sep = "+") %>%
           str_remove_all(fixed(" +")) %>%
           str_remove_all(fixed("+ ")))

(fig_2c_mc <- ggplot(data = control_fluoro_stored, aes(x = wavelength, y = intensity)) +
    geom_line(aes(color = label)) +
    lims(x = c(400, 800)) +
    labs(x = "Wavelength, nm", 
         y = "Intensity, a.u.", 
         color = "Components") +
    theme(legend.text = element_markdown(),
          axis.text.y = element_blank()))

ggsave(plot = fig_2c_mc,
       filename = here("plots", "mendcomm", "fig_2c.tiff"), 
       dpi = 300)




#### order of mixing ####

samples_order <- read_csv(here('data', 'order', 'order.csv')) %>%
  mutate(label = glue("{comp_1_2} < {comp_3} < {comp_4} < {comp_5}")) %>%
  mutate(au_last = ifelse(comp_5 == "Au", TRUE, FALSE)) %>%
  mutate(water_last = ifelse(comp_5 == "w", TRUE, FALSE)) %>%
  mutate(ligands_first = ifelse((comp_4 == "w" & au_last) | comp_4 == "Au" & water_last, TRUE, FALSE))

# import quantum yield fluorescence spectra
order_fluoro <- here('data', 'order', '2023-07-14', 'fluoro') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*-", "", file))) %>%
  left_join(samples_order, by = "file") %>%
  mutate(type = factor(type, levels = c("au", "ref"))) %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  group_by(label, wavelength) %>%
  summarise(intensity = mean(intensity),
            au_last = all(au_last),
            water_last = all(water_last),
            ligands_first = all(ligands_first))



order_abs <- here('data', 'order', '2023-07-14', 'abs') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*-", "", file))) %>%
  left_join(samples_order, by = "file") %>%
  mutate(type = factor(type, levels = c("au", "ref"))) %>%
  filter(absorbance <= 1.9) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  group_by(label, wavelength) %>%
  summarise(absorbance = mean(absorbance),
            au_last = all(au_last),
            water_last = all(water_last),
            ligands_first = all(ligands_first))



# different orders versus all ligands < gold < water
(ggplot(data = order_abs %>% filter(ligands_first), 
       aes(x = wavelength, y = absorbance)) +
  geom_line(aes(group = label), color = "grey", linewidth = 1) +  
  geom_line(data = order_abs %>% filter(ligands_first & water_last),
            aes(group = label), color = "black", linewidth = 1) +
  geom_line(data = order_abs %>% filter(!ligands_first),
            aes(color = label, group = label), linewidth = 1) +
  lims(x = c(282, 500), y = c(0, 1.5)) +
  labs(x = "Wavelength, nm", y = "Absorbance", color = "Mixing order") +
  theme(legend.position = c(0.65, 0.75)))

ggsave(filename = here("plots", "analytica", "fig_1a.tiff"), 
       dpi = 300)


order_abs365 <- order_abs %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(label) %>%
  summarize(absorbance = median(absorbance))

order_fluoro_reduced <- order_fluoro %>%
  left_join(order_abs365, by = "label")



(ggplot(data = order_fluoro_reduced %>% filter(au_last), 
       aes(x = wavelength, y = intensity / absorbance)) +
  geom_line(aes(group = label), color = "black", linewidth = 1) +
  geom_line(data = order_fluoro_reduced %>% filter(water_last),
            aes(group = label), color = "grey", linewidth = 1) +
  geom_line(data = order_fluoro_reduced %>% filter(!ligands_first),
            aes(color = label, group = label), linewidth = 1) +
  lims(x = c(400, 900), y = c(0, 36000)) +
  labs(x = "Wavelength, nm", y = "Normalized intensity", color = "Mixing order") +
    theme(legend.position = c(0.65, 0.75)))

ggsave(filename = here("plots", "analytica", "fig_1b.tiff"), 
       dpi = 300)





#### factorial AuNC ####
# import absorbance spectra for quantum yield measurement
factorial_abs <- here('data', 'factorial', '2023-08-02', 'abs') %>%
  list.files(pattern = "^Absorbance.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type == 'au') %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = glue("au{Au_level}_AMP{AMP_level}_Citr{Citr_level}_NaH{NaH_level}_N3{N3_level}"))



factorial_abs_diluted <- here('data', 'factorial', '2023-08-02', 'diluted_abs') %>%
  list.files(pattern = "^Absorbance.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au', 'ref')) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = glue("au{Au_level}_AMP{AMP_level}_Citr{Citr_level}_NaH{NaH_level}_N3{N3_level}"))

  

# import quantum yield fluorescence spectra
factorial_fluoro_diluted <- here('data', 'factorial', '2023-08-02', 'diluted_fluoro') %>%
  list.files(pattern = "^Subt2_.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au', 'ref')) %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = glue("au{Au_level}_AMP{AMP_level}_Citr{Citr_level}_NaH{NaH_level}_N3{N3_level}"))



factorial_abs365 <- factorial_abs_diluted %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(file) %>%
  summarize(absorbance = median(absorbance))

factorial_qy <- factorial_fluoro_diluted %>%
  left_join(factorial_abs365, by = "file")


# plot for overview of fluorescence spectra
(ggplot(data = factorial_qy %>%
        filter(type == "au") %>%
        filter(wavelength >= 400) %>%
        mutate(intensity = intensity / absorbance)) +
  geom_polygon(data = data.frame(x = c(480, 500, 500, 480),
                                 y = c(-Inf, -Inf, Inf, Inf)),
               aes(x = x, y = y), fill = "turquoise") +
  geom_polygon(data = data.frame(x = c(590, 610, 610, 590),
                                 y = c(-Inf, -Inf, Inf, Inf)),
               aes(x = x, y = y), fill = "gold") +
  geom_line(aes(x = wavelength, y = intensity, group = file, color = factor(N3_level > -1))) +
  labs(x = "Wavelength, nm", y = "Normalzed intensity, a.u.") +
  scale_color_discrete(type = c("blue", "red")) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none"))

ggsave(filename = here("plots", "analytica", "fig_2.tiff"), 
       dpi = 300)


(ggplot(data = factorial_qy %>%
         filter(type == "au") %>%
         filter(N3_level == -1) %>%
         filter(wavelength >= 400) %>%
         mutate(intensity = intensity / absorbance,
                NaH = glue("<i>r</i><sub>Na/H</sub> = {NaH_level}"))) +
  geom_line(aes(x = wavelength, y = intensity, group = file, color = factor(AMP_level))) +
  facet_grid( ~ NaH, scales = "fixed") +
  scale_x_continuous(breaks = seq(400, 800, len = 3)) +
  labs(x = "Wavenumber, nm", y = "Normalized intensity", color = "<i>r</i><sub>AMP</sub> = ") +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "top") +
  theme(strip.text = element_markdown(),
        legend.title = element_markdown()))

ggsave(filename = here("plots", "analytica", "fig_3b.tiff"), 
       dpi = 300)




factorial_fluoro_integral <- factorial_qy %>%
  filter(wavelength > 400) %>%
  group_by(Au_level, AMP_level, Citr_level, NaH_level, N3_level, type, file, label) %>%
  summarize(integral = sum(intensity / absorbance)) %>%
  ungroup()

factorial_scatter_integral <- factorial_qy %>%
  filter(wavelength < 400) %>%
  group_by(Au_level, AMP_level, Citr_level, NaH_level, N3_level, type, file, label) %>%
  summarize(integral = sum(intensity / absorbance)) %>%
  ungroup()



factorial_ref_integral <- factorial_fluoro_integral %>%
  filter(str_detect(type, "ref")) %>%
  select(type, integral) %>%
  rename(ref_integral = integral) %>%
  mutate(type = 'au')


factorial_fluoro_integral <- factorial_fluoro_integral %>%
  left_join(factorial_ref_integral, keep = FALSE, by = "type") %>%
  mutate(qy = integral / ref_integral * 0.54) %>%
  filter(!str_detect(type, "ref"))

factorial_ref_scatter <- factorial_scatter_integral %>%
  filter(str_detect(type, "ref")) %>%
  select(type, integral) %>%
  rename(ref_integral = integral) %>%
  mutate(type = 'au')

factorial_scatter_integral <- factorial_scatter_integral %>%
  left_join(factorial_ref_scatter, keep = FALSE, by = "type") %>%
  mutate(rel_scatter = integral / ref_integral) %>%
  filter(!str_detect(type, "ref"))



factorial_fluoro_intensity_365 <- factorial_qy %>%
  filter(wavelength >= 360 & wavelength <= 370) %>%
  group_by(Au_level, AMP_level, Citr_level, NaH_level, N3_level, type, file, label) %>%
  summarize(int_365 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_365))


factorial_fluoro_intensity_490 <- factorial_qy %>%
  filter(wavelength >= 480 & wavelength <= 500) %>%
  group_by(Au_level, AMP_level, Citr_level, NaH_level, N3_level, type, file, label) %>%
  summarize(int_490 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_490))


factorial_fluoro_intensity_600 <- factorial_qy %>%
  filter(wavelength >= 590 & wavelength <= 610) %>%
  group_by(Au_level, AMP_level, Citr_level, NaH_level, N3_level, type, file, label) %>%
  summarize(int_600 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_600))

factorial_anova <- factorial_fluoro_integral %>%
  select(-c(type, label, integral, ref_integral)) %>%
  left_join(factorial_fluoro_intensity_365, by = "file") %>%
  left_join(factorial_fluoro_intensity_490, by = "file") %>%
  left_join(factorial_fluoro_intensity_600, by = "file") %>%
  mutate(int_600_490 = int_600 / int_490) 





lm_factorial_qy <- lm(data = factorial_anova, 
                  formula = qy * 100 ~ (Au_level + AMP_level + N3_level + Citr_level + NaH_level)^2)
Anova(lm_factorial_qy, type = "II")
summary(lm_factorial_qy)



lm_factorial_qy_reduced <- lm(data = factorial_anova, 
                  formula = qy * 100 ~ N3_level + NaH_level)
Anova(lm_factorial_qy_reduced, type = "II")
summary(lm_factorial_qy_reduced)


# not linear! see plots
lm_factorial_qy_noN3 <- lm(data = factorial_anova %>%
                    filter(N3_level == -1), 
                  formula = qy * 100 ~ Au_level + AMP_level + Citr_level + NaH_level)
Anova(lm_factorial_qy_noN3, type = "II")
summary(lm_factorial_qy_noN3)



lm_factorial_target <- lm(data = factorial_anova %>%
                            filter(qy > 0.01), 
                          formula = int_600_490 ~ (Au_level + AMP_level + N3_level + Citr_level + NaH_level)^2)
Anova(lm_factorial_target, type = "II")
summary(lm_factorial_target)


lm_factorial_target_reduced <- lm(data = factorial_anova %>%
                        filter(qy > 0.01), 
                  formula = int_600_490 ~ Au_level + (AMP_level + N3_level + Citr_level + NaH_level)^2 - AMP_level:Citr_level - AMP_level:NaH_level - Citr_level:NaH_level)
Anova(lm_factorial_target_reduced, type = "II")
summary(lm_factorial_target_reduced)




response_factorial_target <- function(Au, AMP, N3, Citr, NaH){
    coef(lm_factorial_target)["(Intercept)"] + 
    coef(lm_factorial_target)["Au_level"]*Au + 
    coef(lm_factorial_target)["AMP_level"]*AMP + 
    coef(lm_factorial_target)["N3_level"]*N3 + 
    coef(lm_factorial_target)["Citr_level"]*Citr + 
    coef(lm_factorial_target)["NaH_level"]*NaH +
    coef(lm_factorial_target)["AMP_level:N3_level"]*AMP*N3 +
    coef(lm_factorial_target)["N3_level:Citr_level"]*N3*Citr +
    coef(lm_factorial_target)["N3_level:NaH_level"]*N3*NaH}


response_factorial_qy <- function(N3, NaH){
    coef(lm_factorial_qy)["(Intercept)"] + 
    coef(lm_factorial_qy)["N3_level"]*N3 + 
    coef(lm_factorial_qy)["NaH_level"]*NaH}


factorial_predicted <- expand.grid(Au = seq(-1, 1, 0.05), 
                                   AMP = seq(-1, 1, 0.05),
                                   N3 = seq(-1, 1, 0.05),
                                   Citr = seq(-1, 1, 0.05),
                                   NaH = seq(-1, 1, 0.05)) %>%
  mutate(yellow = response_factorial_target(Au = Au, 
                                        AMP = AMP, 
                                        N3 = N3, 
                                        Citr = Citr,
                                        NaH = NaH)) %>%
  mutate(qy = response_factorial_qy(N3 = N3, NaH = NaH)) %>%
  mutate(target = ifelse(qy > 2, qy * yellow, 0))
  


(ggplot(data = factorial_predicted %>%
         filter(Citr %in% c(-1, 0, 1) & NaH %in% c(-1, 0, 1)) %>%
         mutate(NaH = glue("<i>r</i><sub>Na/H</sub> = {NaH}"),
                Citr = glue("<i>r</i><sub>citr</sub> = {Citr}")),
       aes(x = AMP, y = N3, z = qy)) +
  geom_contour_filled() +
  facet_grid(Citr ~ NaH) +
  labs(x = "<i>r</i><sub>AMP</sub>", y = "<i>r</i><sub>azide</sub>", fill = "QY, %") +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme(strip.text.x  = element_markdown(size = 14, face = "bold"),
        strip.text.y  = element_markdown(size = 14, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown()))

ggsave(filename = here("plots", "analytica", "fig_3a.tiff"), 
       dpi = 300)


(ggplot(data = factorial_predicted %>%
         filter(Citr %in% c(-1, 0, 1) & NaH %in% c(-1, 0, 1)) %>%
         mutate(NaH = glue("<i>r</i><sub>Na/H</sub> = {NaH}"),
                Citr = glue("<i>r</i><sub>citr</sub> = {Citr}")),
       aes(x = AMP, y = N3, z = yellow)) +
  geom_contour_filled() +
  facet_grid(Citr ~ NaH) +
  labs(x = "<i>r</i><sub>AMP</sub>", y = "<i>r</i><sub>azide</sub>", 
       fill = "<i>I</i><sub>600</sub>/<i>I</i><sub>490</sub>") +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme(strip.text.x  = element_markdown(size = 14, face = "bold"),
        strip.text.y  = element_markdown(size = 14, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.title = element_markdown()))

ggsave(filename = here("plots", "analytica", "fig_5a.tiff"), 
       dpi = 300)




(ggplot(data = factorial_predicted %>%
         filter(Citr %in% c(-1, 0, 1) & NaH %in% c(-1, 0, 1)) %>%
         mutate(NaH = glue("<i>r</i><sub>Na/H</sub> = {NaH}"),
                Citr = glue("<i>r</i><sub>citr</sub> = {Citr}")),
       aes(x = AMP, y = N3, z = target)) +
  geom_contour_filled() +
  facet_grid(Citr ~ NaH) +
  labs(x = "<i>r</i><sub>AMP</sub>", y = "<i>r</i><sub>azide</sub>", fill = "T, %") +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme(strip.text.x  = element_markdown(size = 14, face = "bold"),
        strip.text.y  = element_markdown(size = 14, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown()))

ggsave(filename = here("plots", "analytica", "fig_5b.tiff"), 
       dpi = 300)





factorial_predicted %>%
  filter(Au == 1) %>%
  arrange(qy) %>%
  tail(100) %>%
  summary()



factorial_predicted %>%
  filter(Au == 1) %>%
  arrange(target) %>%
  tail(20) %>%
  summary()



#### factorial AuNC - additional samples - max yellow ####
# import absorbance spectra for quantum yield measurement
yellow_abs_diluted <- here('data', 'factorial', '2023-08-10', 'diluted_abs') %>%
  list.files(pattern = "^Absorbance.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au3', 'ref')) %>%
  filter(absorbance <= 1.9) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength >= 220) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = glue("au{Au_conc}_AMP{AMP_ratio}_Citr{Citr_ratio}_NaH{NaH_ratio}_N3{N3_ratio}"))


# import quantum yield fluorescence spectra
yellow_fluoro_diluted <- here('data', 'factorial', '2023-08-10', 'diluted_fluoro') %>%
  list.files(pattern = "^Subt2_.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au3', 'ref')) %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = glue("au{Au_conc}_AMP{AMP_ratio}_Citr{Citr_ratio}_NaH{NaH_ratio}_N3{N3_ratio}"))



yellow_abs365 <- yellow_abs_diluted %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(file) %>%
  summarize(absorbance = median(absorbance))

yellow_qy <- yellow_fluoro_diluted %>%
  left_join(yellow_abs365, by = "file")




yellow_fluoro_integral <- yellow_qy %>%
  filter(wavelength > 400) %>%
  group_by(file, Au_conc, AMP_ratio, Citr_ratio, NaH_ratio, N3_ratio, type, label) %>%
  summarize(integral = sum(intensity / absorbance)) %>%
  ungroup() %>%
  filter(!is.na(integral))


yellow_ref_integral <- yellow_fluoro_integral %>%
  filter(str_detect(type, "ref")) %>%
  select(type, integral) %>%
  rename(ref_integral = integral) %>%
  mutate(type = 'au')


yellow_fluoro_integral <- yellow_fluoro_integral %>%
  mutate(ref_integral = mean(yellow_ref_integral$ref_integral)) %>%
  mutate(qy = integral / ref_integral * 0.54) %>%
  filter(!str_detect(type, "ref"))




yellow_fluoro_intensity_365 <- yellow_qy %>%
  filter(wavelength >= 360 & wavelength <= 370) %>%
  group_by(file, Au_conc, AMP_ratio, Citr_ratio, NaH_ratio, N3_ratio, type, label) %>%
  summarize(int_365 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_365))


yellow_fluoro_intensity_490 <- yellow_qy %>%
  filter(wavelength >= 480 & wavelength <= 500) %>%
  group_by(file, Au_conc, AMP_ratio, Citr_ratio, NaH_ratio, N3_ratio, type, label) %>%
  summarize(int_490 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_490))


yellow_fluoro_intensity_600 <- yellow_qy %>%
  filter(wavelength >= 590 & wavelength <= 610) %>%
  group_by(file, Au_conc, AMP_ratio, Citr_ratio, NaH_ratio, N3_ratio, type, label) %>%
  summarize(int_600 = sum(intensity / absorbance / n())) %>%
  ungroup() %>%
  select(c(file, int_600))



yellow_anova <- yellow_fluoro_integral %>%
  select(-c(type, integral, ref_integral)) %>%
  mutate(AMP_level = encode_var(AMP_ratio),
         Citr_level = encode_var(Citr_ratio),
         NaH_level = encode_var(NaH_ratio),
         N3_level = encode_var(N3_ratio)) %>%
  left_join(yellow_fluoro_intensity_365, by = "file") %>%
  left_join(yellow_fluoro_intensity_490, by = "file") %>%
  left_join(yellow_fluoro_intensity_600, by = "file") %>%
  mutate(int_600_490 = int_600 / int_490,
         target =  ifelse(qy > 0.02, qy * int_600_490, 0))




lm_yellow_target <- lm(data = yellow_anova, 
                  formula = 100 * target ~ (AMP_level + N3_level + Citr_level + NaH_level)^2)
Anova(lm_yellow_target, type = "II")
summary(lm_yellow_target)



lm_yellow_target_reduced <- lm(data = yellow_anova, 
                          formula = 100 * target ~ AMP_level + NaH_level + AMP_level:NaH_level + Citr_level)
Anova(lm_yellow_target_reduced, type = "II")
summary(lm_yellow_target_reduced)



lm_yellow_qy <- lm(data = yellow_anova, 
                       formula = qy ~ (AMP_level + N3_level + Citr_level + NaH_level)^2)
Anova(lm_yellow_qy, type = "II")
summary(lm_yellow_qy)


lm_yellow_qy_reduced <- lm(data = yellow_anova, 
                               formula = qy ~ AMP_level + NaH_level)
Anova(lm_yellow_qy_reduced, type = "II")
summary(lm_yellow_qy_reduced)




response_yellow_target = function(AMP, Citr, NaH){
  coef(lm_yellow_target_reduced)["(Intercept)"] + 
    coef(lm_yellow_target_reduced)["AMP_level"]*AMP + 
    coef(lm_yellow_target_reduced)["NaH_level"]*NaH + 
    coef(lm_yellow_target_reduced)["Citr_level"]*Citr + 
    coef(lm_yellow_target_reduced)["AMP_level:NaH_level"]*AMP*NaH}


response_yellow_qy = function(AMP, Citr, NaH){
  coef(lm_yellow_qy_reduced)["(Intercept)"] + 
    coef(lm_yellow_qy_reduced)["AMP_level"]*AMP + 
    coef(lm_yellow_qy_reduced)["NaH_level"]*NaH + 
    coef(lm_yellow_qy_reduced)["Citr_level"]*Citr + 
    coef(lm_yellow_qy_reduced)["AMP_level:NaH_level"]*AMP*NaH}



yellow_predicted <- expand.grid(AMP = seq(-1, 1, 0.05),
                                   Citr = seq(-1, 1, 0.05),
                                   NaH = seq(-1, 1, 0.05)) %>%
  mutate(target = response_yellow_target(AMP = AMP, 
                                     Citr = Citr,
                                     NaH = NaH),
         qy = response_yellow_qy(AMP = AMP, 
                                         Citr = Citr,
                                         NaH = NaH))


(ggplot(data = yellow_predicted %>%
         filter(Citr %in% c(-1, 0, 1)) %>%
         mutate(Citr = glue("<i>r</i><sub>citr</sub> = {Citr}")), 
       aes(x = AMP, y = NaH, z = target)) +
  geom_contour_filled() +
  facet_wrap(~ Citr, nrow = 1) +
  scale_x_continuous(breaks = seq(-1, 1, len = 3)) +
  scale_y_continuous(breaks = seq(-1, 1, len = 3)) +
  labs(x = "<i>r</i><sub>AMP</sub>", y = "<i>r</i><sub>Na/H</sub>", fill = "T, %") +
  theme(strip.text.x  = element_markdown(size = 14, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.position = "top"))

ggsave(filename = here("plots", "analytica", "fig_6.tiff"), 
       dpi = 300)






# yellow_predicted %>%
#   filter(Citr >= 0.95) %>%
#   filter(NaH <= - 0.95) %>%
#   filter(AMP <= -0.95)
# 
# yellow_anova %>%
#   filter(Citr_level >= 0.95) %>%
#   filter(NaH_level <= - 0.95) %>%
#   filter(AMP_level <= -0.95)


#### factorial AuNC - additional samples - max qy ####

# import absorbance spectra for quantum yield measurement
blue_abs_diluted <- here('data', 'factorial', '2023-08-10', 'diluted_abs') %>%
  list.files(pattern = "^Absorbance.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au2', 'ref')) %>%
  filter(absorbance <= 1.9) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength >= 220) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(label = glue("au{Au_conc}_AMP{AMP_ratio}_Citr{Citr_ratio}_NaH{NaH_ratio}_N3{N3_ratio}"))



# import quantum yield fluorescence spectra
blue_fluoro_diluted <- here('data', 'factorial', '2023-08-10', 'diluted_fluoro') %>%
  list.files(pattern = "^Subt2_.*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_factorial, by = "file") %>%
  filter(type %in% c('au2', 'ref')) %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(label = glue("au{Au_conc}_AMP{AMP_ratio}_Citr{Citr_ratio}_NaH{NaH_ratio}_N3{N3_ratio}"))



blue_abs365 <- blue_abs_diluted %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(file) %>%
  summarize(absorbance = median(absorbance))

blue_qy <- blue_fluoro_diluted %>%
  left_join(blue_abs365, by = "file")




blue_fluoro_integral <- blue_qy %>%
  filter(wavelength > 400) %>%
  group_by(file, Au_conc, AMP_ratio, Citr_ratio, NaH_ratio, N3_ratio, type, label) %>%
  summarize(integral = sum(intensity / absorbance)) %>%
  ungroup() %>%
  filter(!is.na(integral))


blue_ref_integral <- blue_fluoro_integral %>%
  filter(str_detect(type, "ref")) %>%
  select(type, integral) %>%
  rename(ref_integral = integral) %>%
  mutate(type = 'au')


blue_fluoro_integral <- blue_fluoro_integral %>%
  mutate(ref_integral = mean(blue_ref_integral$ref_integral)) %>%
  mutate(qy = integral / ref_integral * 0.54) %>%
  filter(!str_detect(type, "ref"))



(ggplot(data = blue_fluoro_integral) +
  geom_col(aes(x = factor(NaH_ratio), y = qy * 100, 
               fill = glue("{AMP_ratio} / {Citr_ratio}")),
           position = position_dodge()) +
  labs(x = "<i>r</i><sub>Na/H</sub>", 
       y = "QY, %", 
       fill = "<i>r</i><sub>AMP</sub> / <i>r</i><sub>citr</sub>") +
  scale_y_continuous(breaks = seq(0, 10, by = 3)) +
  theme(axis.text.x = element_markdown(angle = 45, vjust = 0.7, hjust = 0.5),
        axis.title.x = element_markdown(),
        legend.title = element_markdown(),
        legend.position = "right"))

ggsave(filename = here("plots", "analytica", "fig_4.tiff"), 
       dpi = 300)


#### experiments with samples optimized for qy and target ####
samples_AB_pH <- read_csv(here('data', 'pH', 'samples.csv'))


# import quantum yield fluorescence spectra
AB_pH_fluoro <- here('data', 'pH', '2023-08-24', 'fluoro') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = sub(".*_", "", file)) %>%
  left_join(samples_AB_pH, by = "file") %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(conc_rel = round(V_sample / V_total, 2))



(ggplot(data = AB_pH_fluoro %>% filter(is.na(pH)) %>% filter(conc_rel > 0.1),
        aes(x = wavelength, y = intensity)) +
    geom_line(aes(color = factor(conc_rel))) +
    scale_x_continuous(limits = c(400, 850), breaks = c(400, 600, 800)) +
    labs(x = "Wavelength, nm", 
         y = "Intensity, a.u.", 
         color = "AuNC \nrel. conc.") +
    facet_wrap( ~ sample) +
    theme(legend.title = element_markdown(),
          axis.text.y = element_blank(),
          legend.position = c(0.8, 0.75)))

ggsave(filename = here("plots", "analytica", "fig_7a.tiff"), 
       dpi = 300)



(ggplot(data = AB_pH_fluoro %>% filter(is.na(pH)) %>%
         group_by(file, conc_rel) %>%
         summarize(sample = sample(sample, 1),
                   integral = sum(intensity)),
       aes(x = conc_rel, y = integral, color = sample)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "AuNC rel. conc.", 
       y = "Integral intensity, a.u.", 
       color = "Sample") +
  theme(legend.title = element_markdown(),
        axis.text.y = element_blank(),
        legend.position = c(0.15, 0.8)))

ggsave(filename = here("plots", "analytica", "fig_7b.tiff"), 
       dpi = 300)




AB_pH_abs <- here('data', 'pH', '2023-08-24', 'abs') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = sub(".*_", "", file)) %>%
  left_join(samples_AB_pH, by = "file") %>%
  filter(absorbance <= 1.9) %>%
  filter(wavelength <= 900) %>%
  filter(wavelength <= 650 | wavelength >= 665) %>%
  mutate(conc_rel = round(V_sample / V_total, 2))



AB_pH_abs_365 <- AB_pH_abs %>%
  filter(abs(wavelength - 365) <= 3) %>%
  group_by(file) %>%
  summarize(absorbance = median(absorbance)) %>%
  mutate(absorbance = ifelse(absorbance > 0.01, absorbance, NA))

AB_pH_fluoro_reduced <- AB_pH_fluoro %>%
  left_join(AB_pH_abs_365, by = "file")


  

(ggplot(data = AB_pH_fluoro %>% filter(!is.na(pH)),
        aes(x = wavelength, y = intensity)) +
    geom_line(aes(color = factor(pH))) +
    labs(x = "Wavelength, nm", 
         y = "Intensity, a.u.", 
         color = "pH") +
    scale_x_continuous(limits = c(400, 850), breaks = c(400, 600, 800)) +
    facet_wrap( ~ sample, scales = "free_y") +
    theme(legend.title = element_markdown(),
          axis.text.y = element_blank()))

ggsave(filename = here("plots", "analytica", "fig_8.tiff"), 
       dpi = 300)



# effect of Hg
samples_AB_Hg <- read_csv(here('data', 'Hg', 'samples.csv'))

AB_Hg_fluoro <- here('data', 'Hg', '2023-08-26', 'fluoro') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_AB_Hg, by = "file") %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900) %>%
  mutate(c_Hg = c_Hg_stock * V_Hg / 3300)

(ggplot(data = AB_Hg_fluoro %>%
         filter(wavelength > 400) %>%
          mutate(c_Hg = signif(c_Hg * 1e6, 2)) %>%
          group_by(wavelength, c_Hg, sample, pH) %>%
          summarise(intensity = median(intensity)) %>%
         group_by(c_Hg, sample, pH) %>%
         mutate(intensity = intensity / max(intensity),
                pH = glue("pH = {ifelse(pH == 'low', '3.8', '6.0')}")),
        aes(x = wavelength, y = intensity)) +
    geom_line(aes(color = c_Hg)) +
    labs(x = "Wavelength, nm", 
         y = "Normalised intensity, a.u.", 
         color = "c(Hg), \nμmol/L") +
    scale_x_continuous(limits = c(400, 850), breaks = c(400, 600, 800)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_gradient(low = "red", high = "blue", breaks = c(0, 200, 400)) +
    facet_grid(pH ~ sample) +
    theme(legend.title = element_text(),
          legend.position = "right"))

ggsave(filename = here("plots", "analytica", "fig_10.tiff"), 
       dpi = 300)



AB_Hg_fluoro_integral <- AB_Hg_fluoro %>%
  filter(wavelength > 400) %>%
  group_by(file, c_Hg, sample, pH) %>%
  summarize(integral = sum(intensity)) %>%
  group_by(c_Hg, sample, pH) %>%
  summarize(sd = sd(integral),
            integral = median(integral)) %>%
  ungroup()

AB_Hg_fluoro_ref <- AB_Hg_fluoro_integral %>%
  filter(c_Hg == 0) %>%
  rename(ref = integral) %>%
  select(-c(sd, c_Hg))

AB_Hg_fluoro_integral <- AB_Hg_fluoro_integral %>%
  left_join(AB_Hg_fluoro_ref, by = c("sample", "pH"), keep = FALSE) %>%
  mutate(rel_integral = integral / ref,
         rel_sd = sd / ref)


(ggplot(data = AB_Hg_fluoro_integral %>% filter(c_Hg != 0),
       aes(x = c_Hg, y = rel_integral, color = factor(pH, 
                                                      levels = c("low", "high"), 
                                                      labels = c("3.8", "6.0")))) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE) +
  facet_grid(sample ~ .) +
  scale_x_log10() +
  theme(legend.position = "top") +
  geom_hline(yintercept = 0.97, linetype = 2) +
  labs(x = "c(Hg), mol/L", y = "Rel. intensity", color = "pH") +
  scale_color_discrete(labels = c("3.8", "6.0")))

ggsave(filename = here("plots", "analytica", "fig_9a.tiff"), 
       dpi = 300)



(ggplot(data = AB_Hg_fluoro_integral %>% filter(c_Hg < 1e-4),
       aes(x = c_Hg, y = 1 / rel_integral - 1, color = factor(pH, 
                                                              levels = c("low", "high"), 
                                                              labels = c("3.8", "6.0")))) +
  geom_point() +
  geom_smooth(data = AB_Hg_fluoro_integral %>% filter(c_Hg < 4e-5),
              method = "lm", se = FALSE, formula = y  ~ x - 1, fullrange = TRUE) +
  facet_grid(sample ~ ., scales = "free_y") +
  labs(x = "c(Hg), mol/L", 
       y = expression(paste(frac(italic(I)[0], italic(I)), " \u2212 1")), 
       color = "pH") +
  scale_color_discrete(labels = c("3.8", "6.0")) +
  theme(legend.position = "top",
        axis.title.y = element_text()))

ggsave(filename = here("plots", "analytica", "fig_9b.tiff"), 
       dpi = 300)






# selectivity
samples_AB_selectivity <- read_csv(here('data', 'selectivity', 'samples.csv'))


AB_selectivity_fluoro <- here('data', 'selectivity', '2023-08-28', 'fluoro') %>%
  list.files(pattern = ".*txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro)  %>%
  lapply(FUN = bl_correction, wl_range = c(850, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  mutate(file = as.numeric(sub(".*_", "", file))) %>%
  left_join(samples_AB_selectivity, by = "file") %>%
  filter(wavelength >= 200) %>%
  filter(wavelength <= 900)


AB_Hg_selectivity_integral <- AB_selectivity_fluoro %>%
  filter(wavelength > 400) %>%
  group_by(file, sample, ion) %>%
  summarize(integral = sum(intensity)) %>%
  group_by(sample, ion) %>%
  summarize(sd = sd(integral),
            integral = median(integral)) %>%
  ungroup()

AB_selectivity_fluoro_ref <- AB_Hg_selectivity_integral %>%
  filter(ion == "ref") %>%
  rename(ref = integral) %>%
  select(-c(sd, ion))

AB_Hg_selectivity_integral <- AB_Hg_selectivity_integral %>%
  left_join(AB_selectivity_fluoro_ref, by = "sample", keep = FALSE) %>%
  mutate(rel_integral = integral / ref,
         rel_sd = sd / ref)


(ggplot(data = AB_Hg_selectivity_integral %>% filter(ion != "ref"),
       aes(x = ion, y = rel_integral)) +
  geom_col() +
  facet_grid(sample ~ .) +
  geom_hline(yintercept = 0.97, linetype = 2, color = "red") +
  geom_hline(yintercept = 1.03, linetype = 2, color = "red") +
  theme(axis.text.x = element_markdown(angle = 45, vjust = 0.7, hjust = 0.5)) +
  labs(x = "Probed ion", y = "Rel. intensity"))



ggsave(filename = here("plots", "analytica", "fig_11.tiff"), 
       dpi = 300)










