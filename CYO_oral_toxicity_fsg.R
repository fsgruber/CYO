###########
########### SETUP / Load packages
###########

if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
if(!require(kableExtra)) install.packages("kableExtra", repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) install.packages("gridExtra",repos = "http://cran.us.r-project.org")
if(!require(Rtsne)) install.packages("Rtsne", repos = "http://cran.us.r-project.org")
if(!require(kernlab)) install.packages("kernlab", repos = "http://cran.us.r-project.org")
if(!require(C50)) install.packages("C50",repos = "http://cran.us.r-project.org")
if(!require(ranger)) install.packages("ranger", repos = "http://cran.us.r-project.org")
if(!require(xgboost)) install.packages("xgboost", repos = "http://cran.us.r-project.org")

#limit number of decimal digits displayed
options(digits = 4)


#plotting theme
theme_set(theme_bw())
#personal preference
theme_update(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

#link to the oral toxicity dataset
url <-"https://archive.ics.uci.edu/ml/machine-learning-databases/00508/qsar_oral_toxicity.zip"

#generate temporay file and download data
temp <- tempfile()
download.file(url, temp)

#check what files are in the zip
unzip(temp, list = T)

#extract file and read it in; transform datamode into factor
input_df <- readr::read_delim(unzip(temp), delim = ';', col_names = c(paste0('Bit', seq(1,1024)), 'datamode')) 

#remove temp file
unlink(temp)


###########
########### Data exploration
###########

#check dimensions of input data
input_df %>% dim

#check for na values
input_df %>% is.na() %>% any()

#inspect dataset
input_df %>%
  dplyr::slice(1:10) %>% 
  select(datamode, Bit1:Bit10) %>%
  print()


#visual dataset
set.seed(1998, sample.kind = 'Rounding')
input_df %>% 
  sample_n(50) %>% 
  select(datamode, Bit1:Bit100) %>% 
  mutate(cpd_id = seq(1,50)) %>% 
  gather(key, value, -cpd_id, -datamode) %>% 
  ggplot(aes(key, cpd_id, fill = as.character(value))) +
  geom_tile() +
  geom_point(aes(x = -1, y = cpd_id, color = datamode), size = 0.3) +
  labs(title = 'Oral toxicity data',
       x = 'Fingerprint bit', y = 'Compound',
       subtitle = 'Bit 1: black; Bit 0: white') +
  scale_x_discrete(expand = c(0.1,0)) +
  scale_fill_manual(values = c('white', 'black')) +
  scale_color_manual(name = 'Compound Toxic', values = c('green', 'magenta')) +
  guides(fill = FALSE) +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(color = 'white'),
    plot.subtitle = element_text(face = 'italic'))

#count target classes and visualize
input_df %>% 
  count(datamode) %>% 
  mutate(freq = n/sum(n) * 100,
         freq_lab = paste0(round(freq, digits=0), '%')) %>% 
  ggplot(aes(y = reorder(datamode, n), n)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = freq_lab), hjust = -0.2, size = 3) +
  labs(title = 'Class imbalance', x = 'Number of compounds', y = '') +
  scale_x_continuous(breaks = seq(0,8000,2000), limits = c(0,9100)) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

#perform pca for visualization
pca <- prcomp(input_df %>% select(-datamode))

#visual chemical space
data.frame(pca$x[,1:2], datamode = input_df$datamode) %>% 
  ggplot(aes(PC1, PC2, color = datamode)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(title = 'Chemical Space of dataset') +
  scale_color_manual(name = 'Compound toxic:', values = c('black', 'magenta'))+
  theme(legend.position = 'top')

#calculate percentage variance explained
p1 <- data.frame(variance_explained = pca$sdev^2/sum(pca$sdev^2)*100,
           pc = seq(1, length(pca$sdev))) %>% 
  filter(pc <=50) %>% 
  ggplot(aes(pc, variance_explained)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Principle Component',
       y = '% variance explained')

#plot cumulative sum of percentage variance explained
p2 <- data.frame(variance_explained = pca$sdev^2/sum(pca$sdev^2)*100,
           pc = seq(1, length(pca$sdev))) %>% 
  filter(pc <=50) %>% 
  ggplot(aes(pc, cumsum(variance_explained))) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 55, linetype = 'dashed', color = 'grey') +
  scale_y_continuous(breaks = seq(0,60, 10), limits = c(0,60)) +
  labs(x = 'Principle Component', 
       y = 'Cumulative sum')

#plot both graphs using gridExtra
gridExtra::grid.arrange(p1, p2, ncol = 2)



#perform tsne for visualization
set.seed(1985, sample.kind = 'Rounding')
tsne <- Rtsne(input_df %>% select(-datamode), 
              dims = 3, 
              initial_dims = 50, 
              pca_center = F, 
              check_duplicates = FALSE)

#visualize chemical space using t-sne dimensions
data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], datamode = input_df$datamode) %>% 
  ggplot(aes(tsne1, tsne2, color = datamode)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(name = 'Toxic compound:', values = c('black', 'magenta')) +
  labs(title = expression(~italic(t)~'-SNE visualization of chemical space'),
       x = 't-SNE1', y = 't-SNE2') +
  theme(
    legend.position = 'top'
  )



###########
########### Modeling section
###########

#take only 3 tSNE dimenions for further analysis
new_df <- data.frame(
  tsne$Y
) %>% mutate( y = input_df$datamode %>% factor())


#generate a validation set
set.seed(1998, sample.kind = 'Rounding')
train_index <- createDataPartition(new_df$y, times = 1, p = 0.8, list = FALSE)

#this dataset will be split into training and test set for model evaluation
model_df <- new_df[train_index,]

#this is the hold out set for validating our model
validation_df <- new_df[-train_index,]


#set seed and generate training indices
set.seed(1999, sample.kind = 'Rounding')
train_ind <- createDataPartition(model_df$y, times = 1, p = 0.8, list = FALSE)

train_set <- model_df[train_ind,]
test_set <- model_df[-train_ind,]


### logistic regression model

#ensure downsampling
train_control <- trainControl(sampling = 'down')

#set.seed everytime before training a model to ensure reproducibility
set.seed(1999, sample.kind = 'Rounding')

#train logistic model with all data
fit_glm <- train(y ~ ., data = train_set,
                 method = 'glm',
                 family = 'binomial') 


#train logistic model using downsampling
set.seed(1999, sample.kind = 'Rounding')
fit_glm_ds <- train(y ~ ., data = train_set,
                    method = 'glm',
                    family = 'binomial',
                    trControl = train_control)

#predict compound toxicity
y_hat_glm <- predict(fit_glm, test_set)
y_hat_glm_ds <- predict(fit_glm_ds, test_set)



#performance data frame
perf_df <- tibble(
  Model = 'Logistic Regression all data',
  Accuracy = confusionMatrix(y_hat_glm, test_set$y, positive = 'positive')$overall['Accuracy'],
  Sensitivity = sensitivity(y_hat_glm, test_set$y, positive = 'positive'),
  Specificity = confusionMatrix(y_hat_glm, test_set$y, positive = 'positive')$byClass['Precision'],
  F1 = F_meas(y_hat_glm, test_set$y, relevant = 'positive')
)


#wrapper function to add performance metrics
add_perf <- function(df, model, predictions, true_value) {
  
  df <- df %>% 
    rbind(tibble(
      Model = model,
      Accuracy = confusionMatrix(predictions, true_value, positive = 'positive')$overall['Accuracy'],
      Sensitivity = sensitivity(predictions, true_value, positive = 'positive'),
      Specificity = confusionMatrix(predictions, true_value, positive = 'positive')$byClass['Precision'],
      F1 = F_meas(predictions, true_value, relevant = 'positive')
      
    ))
  
}

#add metric and print
perf_df <- add_perf(perf_df, 'Logistic Regression downsampling', y_hat_glm_ds, test_set$y)
perf_df %>% print()



#### k-nn model 
#set seed, train model
set.seed(1999, sample.kind = "Rounding")
fit_knn <- train(y ~ ., data = train_set,
                method = 'knn',
                trControl = train_control)

#predict compound toxicity
y_hat_knn <- predict(fit_knn, test_set)

#add information to performance data.frame and print 
perf_df <- add_perf(perf_df, 'k-NN', y_hat_knn, test_set$y)
perf_df %>% print()


#### svm model
#set seed and train model
set.seed(1999, sample.kind = 'Rounding')
fit_svm <- train(y ~ ., data = train_set,
                 method = 'svmRadial',
                 trControl = train_control)

#predict compound toxicity
y_hat_svm <- predict(fit_svm, test_set)

#add information to perfomance data.frame and print 
perf_df <- add_perf(perf_df, 'SVM', y_hat_svm, test_set$y)
perf_df %>% print()


#### C5.0 classification tree model
#set seed and train model
set.seed(1999, sample.kind = 'Rounding')
fit_c5 <- train(y ~ ., data = train_set,
                method = 'C5.0',
                trControl = train_control)

#predict compound toxicity
y_hat_c5 <- predict(fit_c5, test_set)

#add information to performance data.frame and print 
perf_df <- add_perf(perf_df, 'C5.0 tree', y_hat_c5, test_set$y)
perf_df %>% print()


##### random forrest model

#set seed
set.seed(1999, sample.kind = 'Rounding')

#use capture.output to suppress warning message; train model
capture.output(
fit_rf <- train(y ~ ., data = train_set,
                method = 'ranger',
                trControl = train_control)
,file ='NUL')

#predict toxicity
y_hat_rf <- predict(fit_rf, test_set)

#add performance metrics to data.frame and print 
perf_df <- add_perf(perf_df, 'random forrest', y_hat_rf, test_set$y)
perf_df %>% print()

#### xgbtree model
#set seed and train model
set.seed(1999, sample.kind = 'Rounding')
fit_xgbtree <- train(y ~ ., data = train_set,
                     method = 'xgbTree',
                     trControl = train_control)

#predict compound toxicity
y_hat_xgbtree <- predict(fit_xgbtree, test_set)

#add performance metrics to data.frame and print
perf_df <- add_perf(perf_df, 'xgb Tree', y_hat_xgbtree, test_set$y)
perf_df %>% print()


#### model ensemble

#create prediction matrix with all model predictions
pred_mat <- tibble(
  
  k_nn = y_hat_knn,
  svm = y_hat_svm,
  rf = y_hat_rf,
  c5 = y_hat_c5,
  xgb = y_hat_xgbtree
  
) %>% as.matrix()


#perform majority vote to predict positive
y_hat_ensemble <- ifelse(rowMeans(pred_mat == 'positive') > 0.5, 'positive', 'negative') %>% factor(levels = levels(test_set$y))

#add performance metrics to data.frame and print
perf_df <- add_perf(perf_df, 'ensemble', y_hat_ensemble, test_set$y)
perf_df %>% print()


###########
########### Model validation
###########

#use random forrest and model ensemble for validation (best sensitivity + f1-score)

#set seed and train random forrest
set.seed(1999, sample.kind = 'Rounding')

#suppress warning message
capture.output(
fit_rf_val <- train(y ~ ., data = model_df,
                    method = 'ranger',
                    trControl = train_control), file = 'NUL')

#predict compounds of validation set
y_hat_rf_val <- predict(fit_rf_val, validation_df)

#add model performance to data.frame and print 
perf_df <- add_perf(perf_df, 'RF validation', y_hat_rf_val, validation_df$y)
perf_df %>% print()



#ensemble validation
models <- c('knn', 'svmRadial', 'C5.0', 'ranger', 'xgbTree')

#suppress warning of RF model
capture.output(
ensemble <- lapply(models, function(model){
  
  set.seed(1999, sample.kind = 'Rounding')
  fit_ens <- train(y ~ ., data = model_df,
                   method = model,
                   trControl = train_control)
  
}), file = 'NUL')

#generate prediction matrix
pred_mat_val <- sapply(ensemble, function(x){
  
  predict(x, validation_df)
  
})

#perform majority vote over models to make ensemble prediction
y_hat_ensemble_val <- ifelse(rowMeans(pred_mat_val == 'positive') > 0.5, 'positive', 'negative') %>% factor(levels(validation_df$y)) 

#add model performance to data.frame and print
perf_df <- add_perf(perf_df, 'ensemble validation', y_hat_ensemble_val, validation_df$y)
perf_df %>% print()

