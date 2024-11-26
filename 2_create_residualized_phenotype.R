setwd("/home/leej5/edwa0506/voting_gwas")
library(tidyverse)
library(readxl)

first_order <- c("vote_2016",
                 "vote_2012",
                 "vote_2008",
                 "vote_2004",
                 "vote_2000",
                 "vote_1996")

second_order <-  c("vote_2018",
                     "vote_2014",
                     "vote_2010",
                     "vote_2006",
                     "vote_2002",
                     "vote_1998",
                     "vote_1994")

all_votes <- c(first_order,
               second_order)

pheno <- read_excel("data/MNVOTE.xlsx")

age_dat <- read_csv("data/ages_x.csv")

pc <- read.table("data/pca.eigenvec")
colnames(pc) <- c("FID", "ID", paste0("PC", 1:20))

data <- left_join(pheno, age_dat, by = "ID") |>
        left_join(pc, by = "ID") |>
        mutate(sex = IDSEX -1, # make it 1 or 0
               age = as.numeric(mdy(BDAY)),
               age2 = age^2,
               age3 = age^3,
               sex_age = sex*age,
               sex_age2 = sex*age2,
               sex_age3 = sex*age3) %>%
        mutate(across(contains("vote"), as.numeric))



all_vote_pheno <- map(1:length(all_votes), \(x) as.formula(paste(all_votes[x]," ~ age + age2 + age3 + sex_age + sex_age2 + sex_age3 +",
                                                        paste("PC", 1:20, sep ="", collapse =" + "), sep=""))) |> 
                  map(\(x) lm(x, data = data, na.action = na.exclude)) |>
                  map(\(x) residuals(x)) |>
                  map(\(x) scale(x)) |>
                  bind_cols() |>
                  rowMeans(na.rm =TRUE) |>
                  na_if(NaN)

second_order_pheno <- map(1:length(second_order), \(x) as.formula(paste(all_votes[x]," ~ age + age2 + age3 + sex_age + sex_age2 + sex_age3 +",
                                                                     paste("PC", 1:20, sep ="", collapse =" + "), sep=""))) |> 
                      map(\(x) lm(x, data = data, na.action = na.exclude)) |>
                      map(\(x) residuals(x)) |>
                      map(\(x) scale(x)) |>
                      bind_cols() |>
                      rowMeans(na.rm =TRUE) |>
                      na_if(NaN)


first_order_pheno <- map(1:length(first_order), \(x) as.formula(paste(all_votes[x]," ~ age + age2 + age3 + sex_age + sex_age2 + sex_age3 +",
                                                                        paste("PC", 1:20, sep ="", collapse =" + "), sep=""))) |> 
                    map(\(x) lm(x, data = data, na.action = na.exclude)) |>
                    map(\(x) residuals(x)) |>
                    map(\(x) scale(x)) |>
                    bind_cols() |>
                    rowMeans(na.rm =TRUE) |>
                    na_if(NaN)

cbind.data.frame(data$IDYRFAM, data$ID, first_order_pheno) |> 
  setNames(c("FID", "IID", "pheno")) |>
  drop_na(pheno) |>
  write.table("data/first_order_pheno.txt", sep = " ", row.names = FALSE, col.names = FALSE)

cbind.data.frame(data$IDYRFAM, data$ID, second_order_pheno) |> 
  setNames(c("FID", "IID", "pheno")) |>
  drop_na(pheno) |>
  write.table("data/first_order_pheno.txt", sep = " ", row.names = FALSE, col.names = FALSE)

cbind.data.frame(data$IDYRFAM, data$ID, all_vote_pheno) |> 
  setNames(c("FID", "IID", "pheno")) |>
  drop_na(pheno) |>
  write.table("data/pheno.txt", sep = " ", row.names = FALSE, col.names = FALSE)

        
