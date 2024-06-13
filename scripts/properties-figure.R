library(ggplot2)
library(dplyr)
library(tidyr)

dat <- read.csv("./properties/calculations_properties_datasets.csv", header = TRUE) %>% 
  setNames(c("id", "dataset", "datatype", "property", "value"))

dataset_names <- unique(dat[["dataset"]])
datatype_names <- unique(dat[["datatype"]])
properties <- unique(dat[["property"]])


cmbns <- combn(dataset_names, 2, simplify = TRUE) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(X1 = factor(X1, levels = unique(X1)),
         X2 = factor(X2, levels = unique(X2)))

type_property <- expand.grid(datatype_names, properties)

all_cmbn_df <- lapply(1L:nrow(type_property), function(ith_row) {
  data.frame(type_property[ith_row, ], cmbns)
}) %>% 
  bind_rows() %>% 
  setNames(c("datatype", "property", "dataset1", "dataset2"))

all_pvals <- lapply(1L:nrow(all_cmbn_df), function(ith_row_id) {
  
  ith_row <- all_cmbn_df[ith_row_id, ]
  
  part_dat <- filter(dat, property == ith_row[["property"]], datatype == ith_row[["datatype"]],
                     dataset %in% unlist(ith_row[, c("dataset1", "dataset2")])) 
  
  res <- try(wilcox.test(value ~ dataset, data = part_dat)[["p.value"]], silent = TRUE)
  
  pval <- if(inherits(res, "try-error")) {
    NA
  } else {
    res
  }
  
  medians <- group_by(part_dat, dataset) %>% 
    summarise(med = median(value, na.rm = TRUE)) %>% 
    pivot_wider(names_from = dataset, values_from = med)
  
  data.frame(pval = pval, 
             median1 = unlist(medians[ith_row[["dataset1"]]], use.names = FALSE), 
             median2 = unlist(medians[ith_row[["dataset2"]]], use.names = FALSE))
}) %>% 
  bind_rows()

plot_df <- cbind(all_cmbn_df, all_pvals) %>% 
  mutate(pval = ifelse(dataset2 == "NP" & datatype != "Full_seq", NA, pval)) %>% 
  group_by(property) %>% 
  mutate(apval = p.adjust(pval, method = "BH", n = sum(!is.na(pval)))) %>% 
  ungroup() %>% 
  mutate(level_signif = cut(apval, breaks = c(0, 0.001, 0.01, 0.05, 1), include.lowest = TRUE, right = TRUE),
         datatype = factor(datatype, levels = c("Full_seq", "IDRs", "PrLDs"),
                           labels = c(expression("Full"~"sequence"), "IDRs", "PrLDs")),
         property = factor(property, levels = c("yr", "kappa_ss", "ncpr", "kappa", "aggrescan", "cars", "anchor"),
                           labels = c(expression("%"~"Y+R"), "kappa[ss]", "NCPR", "kappa", 
                                      expression("Aggregation"~"propensity"),
                                      expression("Cryptic"~"amyloidogenicity"), 
                                      expression("Disorder"~"binding"))),
         level_signif_bigger = factor(ifelse(!is.na(level_signif), paste0(level_signif, 
                                                                          ifelse(median1 > median2, "; 1 > 2", "; 1 < 2")), NA),
                                      levels = c(NA, "[0,0.001]; 1 > 2", "(0.001,0.01]; 1 > 2", "(0.01,0.05]; 1 > 2", "(0.05,1]; 1 > 2", 
                                                 "[0,0.001]; 1 < 2", "(0.001,0.01]; 1 < 2", "(0.01,0.05]; 1 < 2", "(0.05,1]; 1 < 2")))


png("properties-plot.png", width = 1210, height = 640, res = 110)
ggplot(plot_df, aes(x = dataset1, y = dataset2, fill = level_signif)) +
  geom_tile(color = "black") +
  scale_fill_manual("Significance level", values = c(rev(c("#b2e2e2", "#66c2a4", "#238b45")), "#d7301f")) +
  scale_x_discrete("Dataset 1") +
  scale_y_discrete("Dataset 2") +
  facet_grid(datatype ~ property, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

tab_dat <- mutate(plot_df, datatype = factor(datatype, labels = c("Full sequence", "IDRs", "PrLDs")),
                  property = factor(property, labels = c("Y + R%", "kappa ss", "NCPR", "kappa", 
                                                "Aggregation propensity", "Cryptic amyloidogenicity", 
                                                "Disorder binding")),
                  ds = paste0(dataset1, "-", dataset2)) %>% 
  select(ds, datatype, property, apval) %>% 
  mutate(apval = ifelse(apval < 0.01, "X", "")) %>% 
  pivot_wider(id_cols = c(ds, datatype), values_from = apval, names_from = property) %>% 
  select(ds, datatype, `Y + R%`, `kappa ss`, NCPR, kappa, `Aggregation propensity`, 
         `Cryptic amyloidogenicity`, `Disorder binding`)


filter(tab_dat, datatype == "Full sequence", ds %in% c("CE-ND", "DE-ND", "CE-DE")) %>% 
  write.csv(tab_dat, file = "tmp.csv", row.names = FALSE)

