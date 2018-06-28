library(hypervolume) # no working if package rgeos is loaded
library(dplyr)

path <- "~/hypervolume/"
  
environment_trait<- read.csv(file.path(path, "environment_trait.csv"),header = T, dec = ".", sep=",")  
# data brut

environment_trait <- environment_trait %>%  mutate(ID = paste(n_parcelle,n_carre,n_arbre,sep="_"))

environment_trait$LDMC[which(environment_trait$ID == "15_1_198" |environment_trait$ID == "1_2_387")] <- NA

trait_standar<- environment_trait %>% 
  dplyr::select(SLA, LA, LT, LDMC, CC, morphotype_field) %>%
  mutate(morphotype_resample = NA) %>% 
  filter((!is.na(LDMC))) %>% 
  filter(!is.na(LA)) %>% 
  filter(!is.na(LT)) %>%
  filter(!is.na(CC)) %>%
  filter(!is.na(SLA)) # no NA in hypervolume

trait_standar[,1:5] <- as.data.frame(scale(as.matrix(trait_standar[,1:5]),scale = TRUE, center = TRUE))

# tab for result

null_data <- matrix(NA, nrow= 1000, ncol =4)
null_data <- as.data.frame(null_data)
names(null_data) <- c("Method", "G_S","G_SG","S_SG")

# used function local only
.resample <- function(data, ID_name = "morphotype_field", ID_temp = "morphotype_resample"){
  data$morphotype_resample <- data$morphotype_field[sample(1:nrow(data),nrow(data), replace = FALSE)]
  return(data)
}

.calculate_hull <- function(data, ID_morphotype = "morphotype_resample"){
  # subdivise datasets
  data <- data %>% mutate_("morphotype" = ID_morphotype)
  G_standar <-data %>%  filter(morphotype == "G") 
  
  S_standar <-data %>%  filter(morphotype == "S") 
 
  SG_standar <-data %>%  filter(morphotype == "SG") 

  # calculate hull volume
  hull_G <- hypervolume(G_standar[,1:5],method='box')
  hull_S <- hypervolume(S_standar[,1:5],method='box')
  hull_SG <- hypervolume(SG_standar[,1:5],method='box')
  # union and intersection
  
  hull_set1 <- hypervolume_set(hull_G, hull_S, check.memory=FALSE)
  hull_set2 <- hypervolume_set(hull_G, hull_SG, check.memory=FALSE)
  hull_set3 <- hypervolume_set(hull_SG, hull_S, check.memory=FALSE)
  
  # calculate overlapp
  hull_overlap1<- hypervolume_overlap_statistics(hull_set1)
  hull_overlap2<- hypervolume_overlap_statistics(hull_set2)
  hull_overlap3<- hypervolume_overlap_statistics(hull_set3)
  # 
  
  overlap_S_G <- as.data.frame(hull_overlap1)
  overlap_SG_G <- as.data.frame(hull_overlap2)
  overlap_S_SG <- as.data.frame(hull_overlap3)
  
  overlap_total <- cbind(overlap_S_G,overlap_S_SG, overlap_SG_G)
  
  overlap_total <- overlap_total %>% rename(G_S = hull_overlap1, G_SG = hull_overlap2, S_SG = hull_overlap3) %>%
    mutate(Method = row.names(.)) %>% filter(Method == "jaccard")
  
  overlap_total <- overlap_total[, c("Method","G_S", "G_SG", "S_SG")] 
  
  return(overlap_total)
  
}


# rep for 1000 nul obs

for(i in 1:999){
  temp <- .resample(data = trait_standar)
    null_data[i,] <- .calculate_hull(temp)
    print(i)
}

write.csv(null_data, "~/hypervolume/null_data.csv")

