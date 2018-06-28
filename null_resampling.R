null_data <- matrix(NA, nrow= 1000, ncol =4)
null_data <- as.data.frame(null_data)
names(null_data) <- c("Method", "G_S","G_SG","S_SG")
trait_standar <- trait_standar %>% mutate(morphotype_resample = NA)
for(i in 1:999){
  temp <- .resample(data = trait_standar)
    null_data[i,] <- .calculate_hull(temp)
}



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
