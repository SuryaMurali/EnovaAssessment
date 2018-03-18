
# Read in training data set. 
train <- read.csv(file.choose(), sep=',', dec='.', stringsAsFactors = T) 
train$id <- as.character(train$id)
# Changing 'Feb-05' to 2005-02-01
train$diagnosis_date <- as.character(train$diagnosis_date)
mon <- function(x) ifelse(x=='Jan',01,ifelse(x=='Feb',02,ifelse(x=='Mar',03,ifelse(x=='Apr',04,ifelse(x=='May',05,
                   ifelse(x=='Jun',06,ifelse(x=='Jul',07,ifelse(x=='Aug',08,ifelse(x=='Sep',09,ifelse(x=='Oct',10,
                   ifelse(x=='Nov',11,12)))))))))))
train$diagnosis_date <- as.Date(mdy.date(mon(unlist(strsplit(train$diagnosis_date,'-'))[which(c(1:(nrow(train)*2))%%2==1)]),01,as.numeric(paste(20,
                                            unlist(strsplit(train$diagnosis_date,'-'))[which(c(1:(nrow(train)*2))%%2==0)],sep=''))))

#let's see if diagnosis month and year can tell us something
train$diagnosis_month <- as.factor(month(train$diagnosis_date))
train$diagnosis_year <- as.factor(year(train$diagnosis_date))

##from external research: Exposure (uncontrolled) to radiation can cause cell DNA mutations which in turn can cause prostate cancer. Hence adding the variable*
## *Actual data needs to be collected, but for the project's sake, this variable is being randomized.
train$radiation <- 0
set.seed(0827) #for repeatability
train$radiation[sample(1:nrow(train),nrow(train)*0.2,replace = F)] <- 1 #randomly assigning 20% of training data as radiation exposed

for(j in c('race','rd_thrpy','h_thrpy','chm_thrpy','cry_thrpy','brch_thrpy','rad_rem','multi_thrpy','radiation','previous_cancer','smoker')) train[,j] <- as.factor(train[,j])
#'family_history','first_degree_history',
#creating indicator variables for symptoms
train$symptoms <- as.character(train$symptoms)
s <- c()
for(i in 1:nrow(train)) s <- c(s,as.character(unlist(strsplit(train$symptoms[i],','))))
train[,unique(s)] <- 0
for(i in 1:nrow(train)) train[i,as.character(unlist(strsplit(train$symptoms[i],',')))] <- 1
for(j in unique(s)) train[,j] <- as.factor(train[,j])

missmap(train, main = "Missing values vs observed")
## tumor_6_months and psa_6_months have too much missing values, so we remove them.
train$height[which(is.na(train$height))] <- mean(train$height,na.rm=T)
train$weight[which(is.na(train$weight))] <- mean(train$weight,na.rm=T)
ame <- amelia(train,m=5, parallel = "multicore", idvars = c("id",'symptoms','diagnosis_date','tumor_6_months','psa_6_months'),noms=c('t_score', 'n_score', 'm_score', 'stage', 'race',  
                                                                'previous_cancer', 'smoker', 'side', 'rd_thrpy', 
                                                               'h_thrpy', 'chm_thrpy', 'cry_thrpy', 'brch_thrpy', 'rad_rem', 'multi_thrpy', 
                                                               'diagnosis_month', 'diagnosis_year', 'radiation', 'U03', 'U06', 'S07', 'U01', 'U02', 
                                                               'S10', 'O11', 'U05', 'S04', 'P02', 'P01', 'O01', 'O09', 'O08', 'O10', 'P03'))

s <- c()
for(j in names(train)) if(class(train[,j])%in%c('numeric','integer')) {qqnorm(train[,j],main = j)
  qqline(train[,j])}
s <- c()
for(j in names(train)) if(class(train[,j])%in%c('numeric','integer')) s <- c(s,j)
imcdiag(x = train[,s], y =train$survival_1_year)
imcdiag(x = train[,s], y =train$survival_7_years)



##Creating two models (validated by 10 fold Mean Squared Error - measured by Accuracy) 1 yr survival, 7 year survival

#1 year survival
t <- split(c(1:nrow(train)), sample(1:nrow(train), size=10, replace=FALSE))
misClasificError <- vector()
misc <- c()
for (i in 1:10)
{
  for(j in 1:5)
  {
    val <- as.data.frame(ame[[1]][j])
    names(val) <- names(ame[[1]]$imp1)
    val <- val[t[[i]],]
    tr <- as.data.frame(ame[[1]][j])
    names(tr) <- names(ame[[1]]$imp1)
    tr <- tr[-t[[i]],]
    logisticRegression <- glm(survival_1_year ~ gleason_score + t_score + n_score + m_score + stage + age + race + height + weight + family_history + 
                                first_degree_history + previous_cancer + smoker + side + tumor_diagnosis + tumor_1_year + psa_diagnosis +
                                psa_1_year + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                                multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                                U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03
                              , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
    #summary(logisticRegression)
    # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
    # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
    #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
    #Misclassification Error
    fitted.results <- predict(logisticRegression,newdata=val,type='response')
    fitted.results <- ifelse(fitted.results > 0.5,1,0)
    misc[j] <- mean(fitted.results != val$survival_1_year)
  }
  k <- min(which(misc==min(misc)))
  val <- as.data.frame(ame[[1]][k])
  names(val) <- names(ame[[1]]$imp1)
  val <- val[t[[i]],]
  tr <- as.data.frame(ame[[1]][k])
  names(tr) <- names(ame[[1]]$imp1)
  tr <- tr[-t[[i]],]
  logisticRegression <- glm(survival_1_year ~ gleason_score + t_score + n_score + m_score + stage + age + race + height + weight + family_history + 
                              first_degree_history + previous_cancer + smoker + side + tumor_diagnosis + tumor_1_year + psa_diagnosis +
                              psa_1_year + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                              multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                              U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03
                            , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  summary(logisticRegression)
  # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
  # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
  #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #Misclassification Error
  fitted.results <- predict(logisticRegression,newdata=val,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError[i] <- mean(fitted.results != val$survival_1_year)
  misc <- c()
  k <- c()
  # ROC and AREA UNDER CURVE
  p <- predict(logisticRegression,newdata=val,type='response')
  pr <- prediction(p, val$survival_1_year)
  # TPR = sensitivity, FPR=1-specificity
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  plot(prf,main="ROC curve")
  #LIFT CURVE
  perf <- performance(pr,"lift","rpp")
  plot(perf, main="lift curve")
}
misClasificError_1 <- misClasificError

# The MSE is around 10%, so I am proceeding with this model

## final model for 1 year survival involving all training data
misc <- c()
k <- c()
for(j in 1:5)
{
  tr <- as.data.frame(ame[[1]][j])
  names(tr) <- names(ame[[1]]$imp1)
  logisticRegression <- glm(survival_1_year ~ gleason_score + t_score + n_score + m_score + stage + age + race + height + weight + family_history + 
                              first_degree_history + previous_cancer + smoker + side + tumor_diagnosis + tumor_1_year + psa_diagnosis +
                              psa_1_year + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                              multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                              U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03
                            , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #summary(logisticRegression)
  # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
  # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
  #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #Misclassification Error
  fitted.results <- predict(logisticRegression,newdata=val,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misc[j] <- mean(fitted.results != val$survival_1_year)
}
k <- min(which(misc==min(misc)))
tr <- as.data.frame(ame[[1]][k])
names(tr) <- names(ame[[1]]$imp1)
logisticRegression_1 <- glm(survival_1_year ~ gleason_score + t_score + n_score + m_score + stage + age + race + height + weight + family_history + 
                            first_degree_history + previous_cancer + smoker + side + tumor_diagnosis + tumor_1_year + psa_diagnosis + 
                            psa_1_year + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                            multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                            U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03
                          , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
summary(logisticRegression_1)

## survival_7_years

misClasificError <- vector()
misc <- c()
for (i in 1:10)
{
  for(j in 1:5)
  {
    val <- as.data.frame(ame[[1]][j])
    names(val) <- names(ame[[1]]$imp1)
    val <- val[t[[i]],]
    tr <- as.data.frame(ame[[1]][j])
    names(tr) <- names(ame[[1]]$imp1)
    tr <- tr[-t[[i]],]
    
    logisticRegression <- glm(survival_7_years ~ naan(log2(gleason_score)) + t_score + n_score + m_score + stage + naan(log2(age)) + race + naan(log2(height)) + naan(log2(weight)) + family_history + 
                                first_degree_history + previous_cancer + smoker + side + naan(log2(tumor_diagnosis)) +naan(log2(tumor_1_year)) + naan(log2(psa_diagnosis)) + 
                                naan(log2(psa_1_year)) + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                                multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                                U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03 + as.factor(survival_1_year)
                              , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
    #summary(logisticRegression)
    # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
    # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
    #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
    #Misclassification Error
    fitted.results <- predict(logisticRegression,newdata=val,type='response')
    fitted.results <- ifelse(fitted.results > 0.5,1,0)
    misc[j] <- mean(fitted.results != val$survival_7_years)
  }
  k <- min(which(misc==min(misc)))
  val <- as.data.frame(ame[[1]][k])
  names(val) <- names(ame[[1]]$imp1)
  val <- val[t[[i]],]
  tr <- as.data.frame(ame[[1]][k])
  names(tr) <- names(ame[[1]]$imp1)
  tr <- tr[-t[[i]],]
  logisticRegression <- glm(survival_7_years ~ naan(log2(gleason_score)) + t_score + n_score + m_score + stage + naan(log2(age)) + race + naan(log2(height)) + naan(log2(weight)) + family_history + 
                              first_degree_history + previous_cancer + smoker + side + naan(log2(tumor_diagnosis)) +naan(log2(tumor_1_year)) + naan(log2(psa_diagnosis)) + 
                              naan(log2(psa_1_year)) + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                              multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                              U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03 + as.factor(survival_1_year)
                            , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  summary(logisticRegression)
  # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
  # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
  #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #Misclassification Error
  fitted.results <- predict(logisticRegression,newdata=val,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError[i] <- mean(fitted.results != val$survival_7_years)
  misc <- c()
  k <- c()
  # ROC and AREA UNDER CURVE
  p <- predict(logisticRegression,newdata=val,type='response')
  pr <- prediction(p, val$survival_7_years)
  # TPR = sensitivity, FPR=1-specificity
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  plot(prf,main="ROC curve")
  #LIFT CURVE
  perf <- performance(pr,"lift","rpp")
  plot(perf, main="lift curve")
}
misClasificError_7 <- misClasificError

## finam model for survival_7_years
misc <- c()
k <- c()
for(j in 1:5)
{
  tr <- as.data.frame(ame[[1]][j])
  names(tr) <- names(ame[[1]]$imp1)
  logisticRegression <- glm(survival_7_years ~ naan(log2(gleason_score)) + t_score + n_score + m_score + stage + naan(log2(age)) + race + naan(log2(height)) + naan(log2(weight)) + family_history + 
                              first_degree_history + previous_cancer + smoker + side + naan(log2(tumor_diagnosis)) +naan(log2(tumor_1_year)) + naan(log2(psa_diagnosis)) + 
                              naan(log2(psa_1_year)) + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                              multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                              U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03 + as.factor(survival_1_year)
                            , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #summary(logisticRegression)
  # toselect.x <- summary(logisticRegression)$coeff[-1,4] < 0.05
  # logisticRegression <- glm(as.formula(paste('survival_1_year ~ ',paste(names(toselect.x)[toselect.x],collapse=' + '),sep=''))
  #                           , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
  #Misclassification Error
  fitted.results <- predict(logisticRegression,newdata=val,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misc[j] <- mean(fitted.results != val$survival_7_years)
}

k <- min(which(misc==min(misc)))
tr <- as.data.frame(ame[[1]][k])
names(tr) <- names(ame[[1]]$imp1)
logisticRegression_7 <- glm(survival_7_years ~ naan(log2(gleason_score)) + t_score + n_score + m_score + stage + naan(log2(age)) + race + naan(log2(height)) + naan(log2(weight)) + family_history + 
                              first_degree_history + previous_cancer + smoker + side + naan(log2(tumor_diagnosis)) +naan(log2(tumor_1_year)) + naan(log2(psa_diagnosis)) + 
                              naan(log2(psa_1_year)) + tea + rd_thrpy + h_thrpy + chm_thrpy + cry_thrpy + brch_thrpy + rad_rem + 
                              multi_thrpy + diagnosis_month + diagnosis_year + radiation + U03 + U06 + S07 + 
                              U01 + U02 + S10 + O11 + U05 + S04 + P02 + P01 + O01 + O09 + O08 + O10 + P03 + as.factor(survival_1_year)
                            , data = tr, family = binomial(link='logit'),control = list(maxit = 1000))
summary(logisticRegression_7)  

#We have the models now. Let's predict the test data

test <- read.csv(file.choose(),sep=',',dec='.',stringsAsFactors = F)

test$id <- as.character(test$id)
# Changing 'Feb-05' to 2005-02-01
test$diagnosis_date <- as.character(test$diagnosis_date)                                                                                                                                                                                         
test$diagnosis_date <- as.Date(mdy.date(mon(unlist(strsplit(test$diagnosis_date,'-'))[which(c(1:(nrow(test)*2))%%2==1)]),01,as.numeric(paste(20,
                    unlist(strsplit(test$diagnosis_date,'-'))[which(c(1:(nrow(test)*2))%%2==0)],sep=''))))

#let's see if diagnosis month and year can tell us something
test$diagnosis_month <- as.factor(month(test$diagnosis_date))
test$diagnosis_year <- as.factor(year(test$diagnosis_date))

##from external research: Exposure (uncontrolled) to radiation can cause cell DNA mutations which in turn can cause prostate cancer. Hence adding the variable*
## *Actual data needs to be collected, but for the project's sake, this variable is being randomized.
test$radiation <- 0
set.seed(0827) #for repeatability
test$radiation[sample(1:nrow(test),nrow(test)*0.2,replace = F)] <- 1 #randomly assigning 20% of testing data as radiation exposed

for(j in c('race','rd_thrpy','h_thrpy','chm_thrpy','cry_thrpy','brch_thrpy','rad_rem','multi_thrpy','radiation','previous_cancer','smoker')) test[,j] <- as.factor(test[,j])

#creating indicator variables for symptoms
test$symptoms <- as.character(test$symptoms)
s <- c()
for(i in 1:nrow(test)) s <- c(s,as.character(unlist(strsplit(test$symptoms[i],','))))
test[,unique(s)] <- 0
for(i in 1:nrow(test)) test[i,as.character(unlist(strsplit(test$symptoms[i],',')))] <- 1
for(j in unique(s)) test[,j] <- as.factor(test[,j])

missmap(test, main = "Missing values vs observed")

test$height[which(is.na(test$height))] <- mean(test$height,na.rm=T)
test$weight[which(is.na(test$weight))] <- mean(test$weight,na.rm=T)
ame_test <- amelia(test[,which(!(names(test)%in%c('survival_1_year','survival_7_years')))],m=1, parallel = "multicore", idvars = c("id",'symptoms','diagnosis_date','tumor_6_months','psa_6_months'),
                   noms=c('t_score', 'n_score', 'm_score', 'stage', 'race',  
                         'previous_cancer', 'smoker', 'side', 'rd_thrpy', 
                         'h_thrpy', 'chm_thrpy', 'cry_thrpy', 'brch_thrpy', 'rad_rem', 'multi_thrpy', 
                         'diagnosis_month', 'diagnosis_year', 'radiation', 'U03', 'U06', 'S07', 'U01', 'U02', 
                         'S10', 'O11', 'U05', 'S04', 'P02', 'P01', 'O01', 'O09', 'O08', 'O10', 'P03'))

test$survival_1_year_predicted <- predict(logisticRegression_1,newdata=ame_test[[1]]$imp1,type='response')
test$survival_1_year[which(is.na(test$survival_1_year))] <- ifelse(test$survival_1_year_predicted[which(is.na(test$survival_1_year))]>0.5,1,0)
ame_test[[1]]$imp1$survival_1_year <- test$survival_1_year
test$survival_7_years_predicted_prob <- predict(logisticRegression_7,newdata=ame_test[[1]]$imp1,type='response')
test$survival_7_years_predicted <- ifelse(test$survival_7_years_predicted_prob>0.5,1,0)

write.csv(test,'test_scored.csv')
