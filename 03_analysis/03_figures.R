options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

d <- readRDS(file.path("02_processed_data", "pacc_pdac_analysis_cohort.rds"))
fig_dir <- file.path("04_results", "figures")
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(PDAC = "#4C78A8", pACC = "#E6863B")
theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.spacing = grid::unit(7, "pt"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 125, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(paste0(stem, ".pdf"), plot, device = grDevices::cairo_pdf,
         width = width_in, height = height_in, units = "in", bg = "white")
  ggsave(paste0(stem, ".tiff"), plot, device = "tiff", compression = "lzw",
         width = width_in, height = height_in, units = "in", dpi = dpi, bg = "white")
  ggsave(paste0(stem, "_preview.png"), plot, device = "png",
         width = width_in, height = height_in, units = "in", dpi = 180, bg = "white")
  grDevices::svg(paste0(stem, ".svg"), width = width_in, height = height_in,
                 family = "Arial", onefile = TRUE, bg = "white")
  print(plot)
  grDevices::dev.off()
}

# Figure 1: stage-specific OS curves. Regular-grid coordinates preserve the
# step function while keeping editable vector files reasonably compact.
times <- seq(0, 10, by = 0.05)
risk_times <- seq(0, 10, by = 2)
km_rows <- list(); risk_rows <- list(); j <- 0; k <- 0
for (st in c("Localized", "Regional", "Distant")) for (h in c("PDAC", "pACC")) {
  z <- d[d$stage == st & d$histology == h, , drop = FALSE]
  fit <- survfit(Surv(os_time, os_event) ~ 1, data = z, conf.type = "log-log")
  s <- summary(fit, times = times, extend = TRUE)
  j <- j + 1
  km_rows[[j]] <- data.frame(stage = st, histology = h, time = s$time,
                             survival = s$surv, lower = s$lower, upper = s$upper,
                             n_total = nrow(z))
  sr <- summary(fit, times = risk_times, extend = TRUE)
  k <- k + 1
  risk_rows[[k]] <- data.frame(stage = st, histology = h, time = sr$time,
                               n_risk = sr$n.risk)
}
km_plot_data <- do.call(rbind, km_rows)
risk_plot_data <- do.call(rbind, risk_rows)
km_plot_data$stage <- factor(km_plot_data$stage, levels = c("Localized", "Regional", "Distant"))
km_plot_data$histology <- factor(km_plot_data$histology, levels = c("PDAC", "pACC"))
risk_plot_data$stage <- factor(risk_plot_data$stage, levels = c("Localized", "Regional", "Distant"))
risk_plot_data$histology <- factor(risk_plot_data$histology, levels = c("pACC", "PDAC"))
risk_plot_data$hjust <- ifelse(risk_plot_data$time == min(risk_times), 0,
                              ifelse(risk_plot_data$time == max(risk_times), 1, 0.5))
write.csv(km_plot_data, file.path(source_dir, "Figure1_stage_specific_OS.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(risk_plot_data, file.path(source_dir, "Figure1_numbers_at_risk.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p1_curve <- ggplot(km_plot_data, aes(time, survival, colour = histology, fill = histology)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_step(linewidth = 0.65) +
  facet_wrap(~stage, nrow = 1) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = risk_times, limits = c(-0.2, 10.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), labels = function(x) paste0(round(100*x), "%"), expand = c(0, 0)) +
  labs(x = NULL, y = "Overall survival") +
  theme_pub() + theme(legend.position = "top", panel.spacing = grid::unit(12, "pt"), plot.margin = margin(5, 6, 0, 5))

p1_risk <- ggplot(risk_plot_data, aes(time, histology, label = n_risk, colour = histology, hjust = hjust)) +
  geom_text(size = 2.35, show.legend = FALSE) +
  facet_wrap(~stage, nrow = 1) +
  scale_colour_manual(values = palette) +
  scale_x_continuous(breaks = risk_times, limits = c(-0.2, 10.2), expand = c(0, 0)) +
  labs(x = "Time since diagnosis (years)", y = "Number at risk") +
  theme_pub(base_size = 7.5) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    strip.text = element_blank(), legend.position = "none",
    panel.spacing = grid::unit(12, "pt"), plot.margin = margin(0, 6, 5, 5)
  )

p1 <- p1_curve / p1_risk + plot_layout(heights = c(3.5, 1.15))
save_pub(p1, file.path(fig_dir, "Figure1_stage_specific_OS"), height_mm = 125)

# Figure 2: 3-year conditional OS by landmark. Estimates are displayed only
# while at least 20 pACC patients remain at risk in that stage.
cs <- read.csv(file.path("04_results", "tables", "05_conditional_survival.csv"), check.names = FALSE)
cs <- cs[cs$endpoint == "OS" & cs$horizon == 3 & cs$stage %in% c("Localized", "Regional", "Distant"), ]
pacc_risk <- cs[cs$histology == "pACC", c("stage", "landmark", "n_risk")]
names(pacc_risk)[3] <- "pACC_n_risk"
cs <- merge(cs, pacc_risk, by = c("stage", "landmark"), all.x = TRUE)
cs$display <- cs$pACC_n_risk >= 20
cs$stage <- factor(cs$stage, levels = c("Localized", "Regional", "Distant"))
cs$histology <- factor(cs$histology, levels = c("PDAC", "pACC"))
write.csv(cs, file.path(source_dir, "Figure2_conditional_OS.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p2 <- ggplot(cs[cs$display, ], aes(landmark, estimate, colour = histology, group = histology)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.13, linewidth = 0.4) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8, shape = 21, aes(fill = histology), colour = "white", stroke = 0.25) +
  facet_grid(. ~ stage) +
  scale_colour_manual(values = palette) + scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) paste0(round(100*x), "%")) +
  labs(x = "Years already survived", y = "Probability of surviving the next 3 years") +
  guides(fill = "none", colour = guide_legend(override.aes = list(shape = NA, linewidth = 0.8))) +
  theme_pub() + theme(legend.position = "top", strip.text = element_text(colour = "black", face = "bold", margin = margin(b = 4)))
save_pub(p2, file.path(fig_dir, "Figure2_conditional_OS"), height_mm = 82)

# Figure 3: baseline-adjusted histology effect by stage.
cox <- read.csv(file.path("04_results", "tables", "07_cox_histology_effects.csv"), check.names = FALSE)
forest <- cox[cox$model == "Adjusted baseline factors" & cox$endpoint == "OS" &
              cox$comparator == "8140+8500" & cox$stage %in% c("All", "Localized", "Regional", "Distant"), ]
forest$stage_label <- factor(forest$stage, levels = rev(c("All", "Localized", "Regional", "Distant")),
                             labels = rev(c("Overall", "Localized", "Regional", "Distant")))
forest$estimate_label <- sprintf("%.2f (%.2f–%.2f)", forest$HR, forest$lower95, forest$upper95)
write.csv(forest, file.path(source_dir, "Figure3_adjusted_HR.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p3 <- ggplot(forest, aes(HR, stage_label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = 0.18, linewidth = 0.55, colour = palette["pACC"]) +
  geom_point(size = 2.2, shape = 21, fill = palette["pACC"], colour = "white", stroke = 0.3) +
  geom_text(aes(x = 1.12, label = estimate_label), hjust = 0, size = 2.5, family = "sans") +
  scale_x_log10(limits = c(0.30, 1.55), breaks = c(0.4, 0.5, 0.7, 1.0, 1.4)) +
  coord_cartesian(clip = "off") +
  labs(x = "Hazard ratio for pACC vs PDAC (95% CI)", y = NULL) +
  theme_pub() + theme(plot.margin = margin(5, 42, 5, 5))
save_pub(p3, file.path(fig_dir, "Figure3_adjusted_HR"), width_mm = 90, height_mm = 72)

# Figure 4: annual all-cause mortality hazards; late estimates are omitted once
# fewer than 20 pACC patients remain at the interval start.
haz <- read.csv(file.path("04_results", "tables", "06_piecewise_os_hazard.csv"), check.names = FALSE)
haz <- haz[haz$stage %in% c("Localized", "Regional", "Distant"), ]
pacc_haz_risk <- haz[haz$histology == "pACC", c("stage", "interval_start", "n_risk")]
names(pacc_haz_risk)[3] <- "pACC_n_risk"
haz <- merge(haz, pacc_haz_risk, by = c("stage", "interval_start"), all.x = TRUE)
haz$display <- haz$pACC_n_risk >= 20 & haz$hazard_per_100py > 0
haz$midpoint <- (haz$interval_start + haz$interval_end) / 2
haz$stage <- factor(haz$stage, levels = c("Localized", "Regional", "Distant"))
haz$histology <- factor(haz$histology, levels = c("PDAC", "pACC"))
write.csv(haz, file.path(source_dir, "Figure4_piecewise_hazard.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p4 <- ggplot(haz[haz$display, ], aes(midpoint, hazard_per_100py, colour = histology, group = histology)) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.5) +
  facet_wrap(~stage, nrow = 1) +
  scale_colour_manual(values = palette) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
  scale_y_log10() +
  labs(x = "Time since diagnosis (years)", y = "Deaths per 100 person-years (log scale)") +
  theme_pub() + theme(legend.position = "top")
save_pub(p4, file.path(fig_dir, "Figure4_piecewise_hazard"), height_mm = 82)

message("Figures and aggregate source data written to ", fig_dir)
