library(ggplot2)
library(dplyr)
library(tidyr)
library(ranger)
library(mlr3measures)
library(kernlab)

dat <- read.csv("./properties/properties_v5_cars-property.csv", header = FALSE) %>% 
  setNames(c("id", "dataset", "datatype", "property", "value"))

filter(dat, datatype == "Full_seq", property == "kappa") %>% 
  group_by(dataset) %>% 
  summarise(n = length(dataset))

dataset_names <- unique(dat[["dataset"]])
datatype_names <- unique(dat[["datatype"]])
properties <- unique(dat[["property"]])

filter(dat, datatype == "Full_seq") %>% 
  ggplot(aes(x = dataset, y = value)) +
  geom_violin() +
  facet_wrap(~ property, nrow = 1, scales = "free_y")

filter(dat, datatype == "Full_seq") %>% 
  group_by(id, dataset)  %>% 
  summarise(missing_val = sum(is.na(value))) %>% 
  group_by(dataset) %>% 
  summarise(mean(missing_val > 0)) 

non_missing_id <- filter(dat, datatype == "Full_seq") %>% 
  group_by(id, dataset)  %>% 
  summarise(missing_val = sum(is.na(value))) %>% 
  filter(missing_val == 0) %>% 
  pull(id)

only_ful <- filter(dat, datatype == "Full_seq") %>% 
  filter(id %in% non_missing_id)

set.seed(155)

get_folds <- function(x, n) {
  full_folds_abundance <- length(x) %/% n
  full_folds <- rep(1L:n, full_folds_abundance)
  c(sample(full_folds), sample(1L:n, size = length(x) - length(full_folds)))
}

folds <- select(only_ful, id, dataset) %>% 
  unique() %>% 
  group_by(dataset) %>% 
  mutate(fold_id = get_folds(id, 5)) %>% 
  ungroup()

# group_by(folds, dataset, fold_id) %>%
#   summarise(n = length(fold_id)) %>%
#   View()


only_ful_fold <- left_join(only_ful, folds) %>% 
  mutate(dataset = factor(dataset))



train_test_model <- function(train_dat, test_dat, positive) {
  model <- ksvm(x = train_dat[["value"]], y = train_dat[["dataset"]], prob.model = TRUE, type = "C-svc")
  
  preds_df <- data.frame(real = test_dat[["dataset"]],
                         predicted = predict(model, test_dat[["value"]],
                                             type = "probabilities")[, 2]) #%>% 
    #mutate(predicted = factor(predicted, levels = levels(train_dat[["dataset"]])))
  
  reverse_auc(mlr3measures::auc(factor(preds_df[["real"]]), preds_df[["predicted"]], 
                    positive = positive))
}

reverse_auc <- function(x) {
  ifelse(x > 0.5, x, 1 - x)
}

all_folds_res <- lapply(combn(unique(only_ful_fold[["dataset"]]), 2, simplify = FALSE), function(ith_cmb) {
  lapply(1L:5, function(ith_fold) {
    lapply(properties, function(ith_prop) {
      train_dat <- filter(only_ful_fold, property == ith_prop, 
                          datatype == "Full_seq", 
                          fold_id != ith_fold,
                          dataset %in% ith_cmb) %>% 
        droplevels()
      
      test_dat <- filter(only_ful_fold, 
                         property == ith_prop, 
                         datatype == "Full_seq", 
                         fold_id == ith_fold,
                         dataset %in% ith_cmb) %>% 
        droplevels()
      
      data.frame(dataset1 = as.character(ith_cmb[1]),
                 dataset2 = as.character(ith_cmb[2]),
                 fold_id = ith_fold, 
                 property = ith_prop,
                 auc = train_test_model(train_dat, test_dat, positive = as.character(ith_cmb[2])))
    }) %>% bind_rows()
  }) %>% bind_rows()
})  %>% bind_rows() 


all_folds_all_prop_res <- lapply(combn(unique(only_ful_fold[["dataset"]]), 2, simplify = FALSE), function(ith_cmb) {
  lapply(1L:5, function(ith_fold) {
    train_dat <- filter(only_ful_fold, 
                        datatype == "Full_seq", 
                        fold_id != ith_fold,
                        dataset %in% ith_cmb) %>% 
      droplevels() %>% 
      pivot_wider(id_cols = c(id, dataset, datatype, fold_id), names_from = property, values_from = value)
    
    test_dat <- filter(only_ful_fold, 
                       datatype == "Full_seq", 
                       fold_id == ith_fold,
                       dataset %in% ith_cmb) %>% 
      pivot_wider(id_cols = c(id, dataset, datatype, fold_id), names_from = property, values_from = value)
    
    model <- ksvm(dataset ~ kappa + kappa_ss + yr + ncpr + aggrescan + anchor + cars, 
                  data = train_dat, 
                  prob.model = TRUE, type = "C-svc")
    
    preds_df <- data.frame(real = test_dat[["dataset"]],
                           predicted = predict(model, 
                                               test_dat,
                                               type = "probabilities")[,2])  #%>% 
      #mutate(predicted = factor(predicted, levels = levels(train_dat[["dataset"]])))
    
    data.frame(dataset1 = as.character(ith_cmb[1]),
               dataset2 = as.character(ith_cmb[2]),
               fold_id = ith_fold, 
               property = "all properties",
               auc = reverse_auc(mlr3measures::auc(factor(preds_df[["real"]]), preds_df[["predicted"]], 
                                       positive = as.character(ith_cmb[2]))))
  }) %>% bind_rows()
})  %>% bind_rows() 


jpeg("auc.jpeg", width = 1210*1.4, height = 640*1.4, res = 180)
rbind(all_folds_all_prop_res, all_folds_res) %>% 
  group_by(dataset1, dataset2, property) %>% 
  summarise(mean_auc = mean(auc),
            sd_auc = sd(auc)) %>% 
  mutate(property = factor(property, levels = c("yr", "kappa_ss", "ncpr", "kappa", "aggrescan", "cars", "anchor", 
                                                "all properties"),
                           labels = c(expression("%"~"Y+R"), "kappa[ss]", "NCPR", "kappa", 
                                      expression("Aggregation"~"propensity"),
                                      expression("Cryptic"~"amyloidogenicity"), 
                                      expression("Disorder"~"binding"),
                                      expression("All"~"properties")))) %>% 
  mutate(lab = paste0(formatC(mean_auc, digits = 2), "±", formatC(sd_auc, digits = 1))) %>% 
  ggplot(aes(x = dataset1, y = dataset2, fill = mean_auc, label = lab)) +
  geom_tile(color = "black") +
  geom_text(color = "white", size = 2.75) +
  scale_x_discrete("Dataset 1") +
  scale_y_discrete("Dataset 2") +
  scale_fill_gradient("AUC", low = "#df65b0", high = "#67001f") +
  facet_wrap( ~ property, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1))
dev.off()
