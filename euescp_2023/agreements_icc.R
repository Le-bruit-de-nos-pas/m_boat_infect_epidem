

install.packages("vcd")

library("vcd")
library("psych")
library("data.table")
library("tidyverse")


# Categorical Agreement -----------------------------

susceptibilities <- as.table(rbind(
  c(4, 1), c(0,5)
  ))

categories <- c("S", "R")

dimnames(susceptibilities) <- list(Method_1 = categories, Method_2 = categories)

susceptibilities


res.k <- Kappa(susceptibilities)

res.k

confint(res.k)

# Essential Agreement -----------------------------

Method_1 <- data.frame(c(2 , 4 , 8 , 16 , 0.5 , 0.25 , 2 , 4 , 8 , 1))
names(Method_1)[1] <- "Method_1"
Method_2 <- data.frame(c(2 , 16 , 4 , 8, 2 , 1 , 2 , 4 , 4 , 4))
names(Method_2)[1] <- "Method_2"

Methods <- bind_cols(Method_1, Method_2)

ICC(Methods)


# 0.6 C.I. (-0.6 to 0.9)
