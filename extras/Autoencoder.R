#https://statslab.eighty20.co.za/posts/autoencoders_keras_r/

#Consider tying weights https://medium.com/@lmayrandprovencher/building-an-autoencoder-with-tied-weights-in-keras-c4a559c529a2
#https://mc.ai/a-beginners-guide-to-build-stacked-autoencoder-and-tying-weights-with-it/
#Consider 'mae' loss function(mean absolute error) or 'binary_crossentropy'

#' @importFrom zeallot %<-%
#' @import dplyr
#' 
library(dplyr)
epochNum = 500
targetDim = 100
learningRate = 1e-3#1e-2
##preprocess
dim(data) #37689  7395
max(data) #1
min(data) #0
outcomeWeight <- 1/(sum(data)/length(data))
outcomeWeight #35.76255
##Define the encoder and decoder

input_layer <- 
  keras::layer_input(shape = dim(data)[2]) 

lassoRegularizer <- keras::regularizer_l1(l = 1e-4)
ridgeRegularizer <- keras::regularizer_l2(l = 0.01)
elasticRegularizer <- keras::regularizer_l1_l2(l1 = 0.01, l2 = 0.01)

K <- keras::backend()

weighted_mse <- function(y_true, y_pred){
  # convert tensors to R objects
  #y_true   <- K$eval(y_true)
  #y_pred   <- K$eval(y_pred)
  #weights  <- K$eval(outcomeWeight)
  # convert to tensor
  #return(K$constant(loss))
  
  keras::k_mean(keras::k_abs(y_true - y_pred)*keras::k_max(y_true*outcomeWeight, 1))
  
}

weighted_crossentropy <- function(y_true, y_pred){
  keras::k_binary_crossentropy(y_true, y_pred)*keras::k_max( y_true*outcomeWeight, 1)
  #keras::k_sum(-((y_true*keras::k_log(keras::k_abs(y_pred)))+(1-y_true)*keras::k_log(keras::k_abs(1-y_pred))))*keras::k_max(y_true*outcomeWeight, 1)
}

metric_weighted_mse <- keras::custom_metric("weighted_mse", function(y_true, y_pred) {
  weighted_mse(y_true, y_pred)
})

metric_f1 <- function (y_true,y_pred) {
  y_pred <- keras::k_round(y_pred)
  precision <- keras::k_sum(y_pred*y_true)/(keras::k_sum(y_pred)+keras::k_epsilon())
  recall    <- keras::k_sum(y_pred*y_true)/(keras::k_sum(y_true)+keras::k_epsilon())
  (2*precision*recall)/(precision+recall+keras::k_epsilon())
} 

# encoder <-
#   input_layer %>%
#    keras::layer_dense(units = targetDim*2, activation = "relu", activity_regularizer = lassoRegularizer
#                       ) %>%
#   keras::layer_batch_normalization() %>%
#   keras::layer_dropout(rate = 0.2) %>%
  # keras::layer_dense(units = 1000, activation = "relu", kernel_regularizer = NULL) %>%
  # keras::layer_dropout(rate = 0.1) %>%
  # keras::layer_dense(units = 250, activation = "relu", kernel_regularizer = NULL) %>%
#  keras::layer_dense(units = targetDim, activity_regularizer = NULL,activation = "sigmoid") # 2 dimensions for the output layer

encoder <-
  input_layer %>%
  keras::layer_dense(units = targetDim, activity_regularizer = lassoRegularizer,activation = "sigmoid") # 2 dimensions for the output layer

##Consider weight regularization 
#

decoder <- 
  encoder %>% 
  # keras::layer_dense(units = targetDim*2, activation = "relu", kernel_regularizer = NULL) %>% 
  # keras::layer_dropout(rate = 0.2) %>% 
  # keras::layer_dense(units = 1000, activation = "relu", kernel_regularizer = NULL) %>%
  # keras::layer_dropout(rate = 0.1) %>%
  # keras::layer_dense(units = 5000, activation = "sigmoid", kernel_regularizer = NULL) %>%
  keras::layer_dense(units = dim(data)[2], activation = "relu") #the original dimension

##compile and train the autoencoder
autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

autoencoder_model %>% keras::compile(
  loss=weighted_crossentropy,#'mean_squared_error',#weighted_mse
  optimizer= keras::optimizer_adam(lr = learningRate),#lr=1e-2),
  metrics = c("binary_crossentropy", "binary_accuracy")
)

earlyStopping = keras::callback_early_stopping(monitor = "val_loss", 
                                               patience = 40,
                                               mode="auto",
                                               min_delta = 0)
reduceLr = keras::callback_reduce_lr_on_plateau(monitor="val_loss", factor =0.1, 
                                              patience = 15,mode = "auto", min_delta = 1e-5, cooldown = 0, min_lr = 0)
#summary(autoencoder_model)


##train onto itself
history <-
  autoencoder_model %>%
  keras::fit(data,
             data,
             epochs = epochNum,
             shuffle = TRUE,
             validation_split = 0.1,
             callbacks = list(earlyStopping, reduceLr)
             )

autoencoder_model %>% keras::save_model_hdf5(file.path(exportFolder,"autoencoder_model.h5"))
# autoencoder_weights <- 
#   autoencoder_model %>%
#   keras::get_weights()
autoencoder_model %>% keras::save_model_weights_hdf5(file.path(exportFolder,"autoencoder_model_weights.h5"),overwrite = TRUE)


####Load the weights####
# autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)
# autoencoder_model %>% keras::load_model_weights_hdf5(file.path(exportFolder,"autoencoder_model_weights.h5"), skip_mismatch = TRUE, by_name = TRUE)
# autoencoder_model %>% keras::compile(
#   loss=weighted_crossentropy,#'mean_squared_error',#weighted_mse
#   optimizer= keras::optimizer_adam(),#lr=1e-2),
#   metrics = c("binary_crossentropy", "binary_accuracy")
# )
# 
encoder_model <- keras::keras_model(inputs = input_layer, outputs = encoder)
encoder_model %>% keras::load_model_weights_hdf5(file.path(exportFolder,"autoencoder_model_weights.h5"), skip_mismatch = TRUE, by_name = TRUE)
encoder_model %>% keras::compile(
  loss= weighted_crossentropy,#'mean_squared_error',#weighted_mse,#'mean_squared_error',#weighted_mse
  optimizer= keras::optimizer_adam(lr=learningRate),
  metrics = c("binary_crossentropy", "accuracy",metric_f1)
)


#0.6135
#binary_crossentropy: 1.4176 - binary_accuracy: 0.9058 - val_loss: 31.9301 - val_binary_crossentropy: 1.4257 -val_binary_accuracy: 0.9053
#loss: 15.0099 - binary_crossentropy: 1.3001 - binary_accuracy: 0.9135 - val_loss: 15.5911 - val_binary_crossentropy: 1.3106 - val_binary_accuracy: 0.9130
#loss: 2.6108 - binary_crossentropy: 0.0785 - binary_accuracy: 0.9756 - val_loss: 2.8138 - val_binary_crossentropy: 0.0809 - val_binary_accuracy: 0.9748
#loss: 2.6155 - binary_crossentropy: 0.0787 - binary_accuracy: 0.9763 - val_loss: 2.8157 - val_binary_crossentropy: 0.0809 - val_binary_accuracy: 0.9755
#loss: 5.6514 - binary_crossentropy: 0.0793 - binary_accuracy: 0.9917 - val_loss: 6.4818 - val_binary_crossentropy: 0.0909 - val_binary_accuracy: 0.9904
#Epoch200 loss: 6.2063 - binary_crossentropy: 0.1765 - binary_accuracy: 0.9841 - val_loss: 6.2436 - val_binary_crossentropy: 0.1787 - val_binary_accuracy: 0.9838
#Epoch200 loss: 6.0014 - binary_crossentropy: 0.1690 - binary_accuracy: 0.9876 - val_loss: 6.0218 - val_binary_crossentropy: 0.1706 - val_binary_accuracy: 0.9870

#autoencoder_model <- keras::load_model_hdf5(file.path(exportFolder,"autoencoder_model.h5"))
#encoder_model <- keras::load_model_hdf5(file.path(exportFolder,"encoder_weight.h5"))


encodedData <- autoencoder_model %>% 
  keras::predict_on_batch (data[1:10000,])

data[1:10,20:30]
encodedData[1:10,20:30]

reducedDim <- encoder_model %>% 
  keras::predict_on_batch (data)

#MSE
mean((as.matrix(data[1:10000,]) -encodedData)^2) # 0.00621386
#Absolute difference mean
mean(abs(as.matrix(data[1:10000,]) -encodedData)) # 0.01516118

#MSE of non-zero values
index <- which(as.matrix(data[1:10000,])==1)
mean((as.matrix(data[1:10000,])[index] -encodedData[index])^2) # 0.3452044
#Absolute difference mean  of non-zero values
mean(abs(as.matrix(data[1:10000,])[index] -encodedData[index])) # 0.4550976


#target dim = 100 / 200 / 100(epoch 500)

####AUROC
pred <- ROCR::prediction(c(encodedData), c(as.matrix(data[1:10000,])))
ROCR::performance(pred,"auc")@y.values[[1]] # 0.8046417 / 0.8065142 /  0.8190414

#perfPrc <- ROCR::performance(pred,"prec","rec")

####Sensitivity & Specificity
perfRoc <- ROCR::performance(pred,"tpr","fpr")
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
  return(cut.ind)
}
(opt.cut(perfRoc, pred)*100)[1] #sensitivity #62.01026 / 61.57023 / 64.55205
(opt.cut(perfRoc, pred)*100)[2] #specificity #97.64553 / 99.3522 / 98.27794

####AUPRC
prcCurv <- PRROC::pr.curve(scores.class0 = c(encodedData), weights.class0=c(as.matrix(data[1:10000,])),curve=TRUE)
#plot(prcCurv)
prcCurv$auc.integral #0.6036849 / 0.638942 /  0.6455071
