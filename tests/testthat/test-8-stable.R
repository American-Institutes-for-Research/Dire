if(FALSE) {
source("stable_helpers.R")
#this will be fixed, item par, block design 





set.seed(345)
n_items <- 5
N <- 10000
#correlation matrix 6 variables
d <- matrix(data = .4, nrow = 6, ncol = 6)
diag(d) <- 1
c_prop <- list(c(1),
               c(.4,1),
               c(.1,.3,.5,.7,1),
               c(.3,.6,1),
               c(.7,.9,1),
               c(.4,.6,.9,1))
#30 items, these are item charactheristics
item <- as.integer(c(1:n_items))
a <- runif(n_items, 0.5, 1.5)
b <- runif(n_items, -3, 3)
c <- runif(n_items, 0, .15)
item_par <- data.frame(item, a, b, c)
#no block design
block_ex <- block_design(n_blocks = 1, item_parameters = item_par)
blk_book <- matrix(1, nrow = 1, ncol = 1)
book_ex <- booklet_design(item_block_assignment = block_ex$block_assignment, book_design = blk_book)
book_admin <- booklet_sample(n_subj = N, book_item_design = book_ex)
#5 background var and 1 theta
bgr_dt <- questionnaire_gen(n_obs = N, cat_prop = c_prop, cor_matrix = d, n_vars = 6, theta = TRUE) 
bgr_dt$q1 <- ifelse(bgr_dt$q1 ==1, 0 ,1) #turn into zer0 one instead of 1-2
resp_matrix <- response_gen2(subject = book_admin$subject, item = book_admin$item, theta = bgr_dt$theta, a_par = a, b_par = b, c_par = c)

names(item_par) <- c("ItemID", "slope", "difficulty", "guessing")
item_par$test <- "test1"
item_par$subtest <- "subtest1"
item_par$D <- 1
item_par$MODEL <- "3pl"
item_names <- names(resp_matrix[,1:n_items])
item_par$ItemID <- item_names
dichotParamTab <- item_par[,c("ItemID", "test", "subtest", "slope", "difficulty", "guessing", "D", "MODEL" )]
#str(resp_matrix)
subject <- resp_matrix$subject #vector of student ids
stuItems <- reshape(data=resp_matrix, varying=c(item_names), idvar=c("subject"),
        direction="long", v.names="score", times=item_names, timevar="key")
stuDat_baseline <- data.frame(subject, q1=bgr_dt$q1)#this will change based on the model
mml_baseline <- mml(test1 ~ q1, stuItems = stuItems, stuDat = stuDat_baseline, dichotParamTab=dichotParamTab, idVar="subject", Q=34)
summary(mml_baseline)
summary(lm(theta ~ q1, data=bgr_dt))
# agrees







set.seed(345)
n_items <- 30
N <- 10000
#correlation matrix 6 variables
d <- matrix(data = .4, nrow = 6, ncol = 6)
diag(d) <- 1
c_prop <- list(c(1),
               c(.4,1),
               c(.1,.3,.5,.7,1),
               c(.3,.6,1),
               c(.7,.9,1),
               c(.4,.6,.9,1))
#30 items, these are item charactheristics
item <- as.integer(c(1:n_items))
a <- runif(n_items, 0.5, 1.5)
b <- runif(n_items, -3, 3)
c <- runif(n_items, 0.05,0.2)
item_par <- data.frame(item, a, b, c)
#no block design
block_ex <- block_design(n_blocks = 1, item_parameters = item_par)
blk_book <- matrix(1, nrow = 1, ncol = 1)
book_ex <- booklet_design(item_block_assignment = block_ex$block_assignment, book_design = blk_book)
book_admin <- booklet_sample(n_subj = N, book_item_design = book_ex)
#5 background var and 1 theta
bgr_dt <- questionnaire_gen(n_obs = N, cat_prop = c_prop, cor_matrix = d, n_vars = 6, theta = TRUE) 
bgr_dt$q1 <- ifelse(bgr_dt$q1 ==1, 0 ,1) #turn into zer0 one instead of 1-2
resp_matrix <- response_gen2(subject = book_admin$subject, item = book_admin$item, theta = bgr_dt$theta, a_par = a, b_par = b, c_par = c)

bb <- cbind(bgr_dt, resp_matrix)
bb <- bb[order(bb$theta),]
bb$thetai <- cut(bb$theta,breaks=c(-Inf, quantile(bb$theta,seq(0.01,0.99,len=99)),Inf))
for(j in unique(bb$thetai)) {
  bb$thetaj[bb$thetai==j] <- mean(bb$theta[bb$thetai==j])
}
for(i in 1:ncol(resp_matrix)) {
  iname <- colnames(resp_matrix)[i]
  plot(range(bb$thetaj), c(0,1), type="n")
  bb$thisq <- bb[,iname]
  minidat <- aggregate(thisq ~ thetaj, data=bb, FUN=mean)
  points(minidat$thetaj, minidat$thisq, xlab="theta", ylab="frac correct")
  yt <- c[i] + (1-c[i]) /(1+exp(-1*a[i]*(minidat$thetaj - b[i])))
  lines(minidat$thetaj, yt, col="red")
#  Sys.sleep(5)
}


names(item_par) <- c("ItemID", "slope", "difficulty", "guessing")
item_par$test <- "test1"
item_par$subtest <- "subtest1"
item_par$D <- 1
item_par$MODEL <- "3pl"
item_names <- names(resp_matrix[,1:n_items])
item_par$ItemID <- item_names
dichotParamTab <- item_par[,c("ItemID", "test", "subtest", "slope", "difficulty", "guessing", "D", "MODEL" )]
#str(resp_matrix)
subject <- resp_matrix$subject #vector of student ids
stuItems <- reshape(data=resp_matrix, varying=c(item_names), idvar=c("subject"),
        direction="long", v.names="score", times=item_names, timevar="key")
stuDat_baseline <- data.frame(subject, q1=bgr_dt$q1)#this will change based on the model

mml_baseline <- mml(test1 ~ q1, stuItems = stuItems, stuDat = stuDat_baseline, dichotParamTab=dichotParamTab, idVar="subject", Q=34)
summary(mml_baseline)
summary(lm(theta ~ q1, data=bgr_dt))
# agrees



n_items <- 99
N <- 10000
#correlation matrix 6 variables
d <- matrix(data = .4, nrow = 6, ncol = 6)
diag(d) <- 1
c_prop <- list(c(1),
               c(.4,1),
               c(.1,.3,.5,.7,1),
               c(.3,.6,1),
               c(.7,.9,1),
               c(.4,.6,.9,1))
#30 items, these are item charactheristics
item <- as.integer(c(1:n_items))
a <- runif(n_items, 0.5, 1.5)
b <- runif(n_items, -3, 3)
c <- runif(n_items, 0, .15)
item_par <- data.frame(item, a, b, c)
#no block design
block_ex <- block_design(n_blocks = 1, item_parameters = item_par)
blk_book <- matrix(1, nrow = 1, ncol = 1)
book_ex <- booklet_design(item_block_assignment = block_ex$block_assignment, book_design = blk_book)
book_admin <- booklet_sample(n_subj = N, book_item_design = book_ex)
#5 background var and 1 theta
bgr_dt <- questionnaire_gen(n_obs = N, cat_prop = c_prop, cor_matrix = d, n_vars = 6, theta = TRUE) 
bgr_dt$q1 <- ifelse(bgr_dt$q1 ==1, 0 ,1) #turn into zer0 one instead of 1-2
resp_matrix <- response_gen2(subject = book_admin$subject, item = book_admin$item, theta = bgr_dt$theta, a_par = a, b_par = b, c_par = c)

names(item_par) <- c("ItemID", "slope", "difficulty", "guessing")
item_par$test <- "test1"
item_par$subtest <- "subtest1"
item_par$D <- 1
item_par$MODEL <- "3pl"
item_names <- names(resp_matrix[,1:n_items])
item_par$ItemID <- item_names
dichotParamTab <- item_par[,c("ItemID", "test", "subtest", "slope", "difficulty", "guessing", "D", "MODEL" )]
#str(resp_matrix)
subject <- resp_matrix$subject #vector of student ids
stuItems <- reshape(data=resp_matrix, varying=c(item_names), idvar=c("subject"),
        direction="long", v.names="score", times=item_names, timevar="key")
stuDat_baseline <- data.frame(subject, q1=bgr_dt$q1)#this will change based on the model
mml_baseline <- mml(test1 ~ q1, stuItems = stuItems, stuDat = stuDat_baseline, dichotParamTab=dichotParamTab, idVar="subject", Q=34)
summary(mml_baseline)
summary(lm(theta ~ q1, data=bgr_dt))


n_items <- 500
N <- 10000
#correlation matrix 6 variables
d <- matrix(data = .4, nrow = 6, ncol = 6)
diag(d) <- 1
c_prop <- list(c(1),
               c(.4,1),
               c(.1,.3,.5,.7,1),
               c(.3,.6,1),
               c(.7,.9,1),
               c(.4,.6,.9,1))
#30 items, these are item charactheristics
item <- as.integer(c(1:n_items))
a <- runif(n_items, 0.5, 1.5)
b <- runif(n_items, -3, 3)
c <- runif(n_items, 0, .15)
item_par <- data.frame(item, a, b, c)
#no block design
block_ex <- block_design(n_blocks = 1, item_parameters = item_par)
blk_book <- matrix(1, nrow = 1, ncol = 1)
book_ex <- booklet_design(item_block_assignment = block_ex$block_assignment, book_design = blk_book)
book_admin <- booklet_sample(n_subj = N, book_item_design = book_ex)
#5 background var and 1 theta
bgr_dt <- questionnaire_gen(n_obs = N, cat_prop = c_prop, cor_matrix = d, n_vars = 6, theta = TRUE) 
bgr_dt$q1 <- ifelse(bgr_dt$q1 ==1, 0 ,1) #turn into zer0 one instead of 1-2
resp_matrix <- response_gen2(subject = book_admin$subject, item = book_admin$item, theta = bgr_dt$theta, a_par = a, b_par = b, c_par = c)

names(item_par) <- c("ItemID", "slope", "difficulty", "guessing")
item_par$test <- "test1"
item_par$subtest <- "subtest1"
item_par$D <- 1
item_par$MODEL <- "3pl"
item_names <- names(resp_matrix[,1:n_items])
item_par$ItemID <- item_names
dichotParamTab <- item_par[,c("ItemID", "test", "subtest", "slope", "difficulty", "guessing", "D", "MODEL" )]
#str(resp_matrix)
subject <- resp_matrix$subject #vector of student ids
stuItems <- reshape(data=resp_matrix, varying=c(item_names), idvar=c("subject"),
        direction="long", v.names="score", times=item_names, timevar="key")
stuDat_baseline <- data.frame(subject, q1=bgr_dt$q1)#this will change based on the model
mml_baseline <- mml(test1 ~ q1, stuItems = stuItems, stuDat = stuDat_baseline, dichotParamTab=dichotParamTab, idVar="subject", Q=34)
summary(mml_baseline)
summary(lm(theta ~ q1, data=bgr_dt))



}