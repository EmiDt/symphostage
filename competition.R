# Packages devtools et pingr nécessaires
# Package EcoFoG pour l'accès aux données de Paracou
# devtools::install_github("EcoFoG/EcoFoG")
library("EcoFoG")
# Package dbmss (extension de spatstat)
library("dbmss")
library("tidyverse")
library("pingr")
library("spatstat")
library("devtools")
library("SpatDiv")


# chargement des données manuelle 


paracoutree_raw <- read.csv("C:\\Users\\emduc\\Desktop\\Drive\\symphostage\\DataParacou_geraldine.csv",
                        header = T, dec = ".", sep=",")

symphonia_stage <- read.csv("C:\\Users\\emduc\\Desktop\\Drive\\symphostage\\data_merging\\Full_traits_data.csv",
                            header = T, dec = ",", sep=";")%>% 
  group_by(n_parcelle, n_carre, n_arbre) %>% 
  summarize() %>% 
  mutate(ID = paste(n_parcelle, n_carre, n_arbre, sep="_"))

paracoutree <- paracoutree_raw %>%
  filter( Plot %in% 1:16) %>% 
  filter(CensusYear == "2015") %>% 
  dplyr::select(Plot, SubPlot, TreeFieldNum, Xfield, Yfield, Genus, Species,
                CircCorr, CodeAlive ,Xutm, Yutm ) %>% 
  rename(n_parcelle ="Plot", n_carre = "SubPlot", n_arbre ="TreeFieldNum",
  X = "Xfield", Y="Yfield", circ_corr = "CircCorr", code_vivant = "CodeAlive") %>%
  mutate(species = paste(Genus,Species, sep="_")) 





# lBAL --------------------------------------------------------------------


indices = matrix(ncol = 4)
indices <- as.data.frame(indices); names(indices) <- c("n_parcelle","n_carre", "n_arbre","GNeighbors")


for (i in 1:15) {
  
Paracoutest <- paracoutree %>% 
  filter(code_vivant ==1) %>% 
  select(n_parcelle, n_carre, n_arbre, X, Y, circ_corr, species) %>% 
  filter(n_parcelle == i) %>% 
  mutate( PointWeight = (circ_corr)^2/(4*pi)) %>% 
  mutate( PointName = paste(n_parcelle, n_carre, n_arbre, sep="_")) %>%
  rename( PointType = species) %>% 
  na.omit()
rownames(Paracoutest) <- Paracoutest$PointName
Paracoutest <- dbmss::wmppp(Paracoutest, window = owin(c(0,250), c(0,250),
                             unitname=c("metres", "metres")))

# Définition du voisinage

rVoisinage <- 25


# 1 : sans correction des effets de bord
###########
# Matrice de distance
Distances <- pairdist(Paracoutest)
# Matrice de voisinage (25m)
Voisins <- (Distances <= rVoisinage) # inférieur ou égale


  # # Variante : seulement les gros arbres
  # Voisins <- (Distances <= rVoisinage) & (Paracoutest$marks$PointWeight > 1000)
  # 
  # # Variante 2 : les symphonia
  # EstSymphonia <- logical(Paracoutest$n) #vecteur de la longueur de paracoutest ici 3540
  # EstSymphonia[grep("Symphonia_", Paracoutest$marks$PointType)] <- TRUE #Èremplacement par TRUE des place de dans le vecteur qui
  # #corresponde à celle dans lesquelles il y a le pattern Symphonia_
  # Voisins <- (Distances <= rVoisinage) & EstSymphonia # on selectionne dans la matrice voisin ceux qui répondent aux 2 conditions
  
  
# Elimination du point lui-même
diag(Voisins) <- FALSE

# Surface terrière des voisins, par colonne

Gvoisins <- apply(Voisins, 2, function(EstVoisin) sum(Paracoutest[EstVoisin]$marks$PointWeight))


# 2 : correction des effets de bord
###########

# Facteur de correction
Correction <- function(NumPoint) {
  # Disque de 25m de rayon autour du point
  disc(radius=rVoisinage, centre=c(Paracoutest$x[NumPoint], Paracoutest$y[NumPoint])) %>% 
  # Intersection avec la parcelle
  intersect.owin(Paracoutest$window) %>% 
  # Calul de la surface
  area -> VoisinageDansParcelle
  # Retour du facteur de correction
  return(pi * rVoisinage^2 / VoisinageDansParcelle)
}

Corrections <- vapply(1:Paracoutest$n, Correction, 0)

# Surface terrière corrigée ajoutée au wmppp
Paracoutest$marks$GNeighbors <- Gvoisins*Corrections


# Nombre de voisins
#Paracoutest$marks$nNeighbors <- colSums(Voisins)*Corrections


# # Carte
# plot(density(Paracoutest, weights = Paracoutest$marks$GNeighbors), main="Surface terrière des voisins")
# plot(Paracoutest[grep("Symphonia_", Paracoutest$marks$PointType)], which.marks = "PointWeight", add=TRUE)
# 
# # Surface terrière
# plot(density(Paracoutest, weights = Paracoutest$marks$PointWeight), main="Surface terrière")

# extraction des données 

P_Gneighbour <- as.data.frame(Paracoutest) 
P_Gneighbour <- P_Gneighbour %>% 
  mutate(ID = rownames(P_Gneighbour)) 
 G_neighbour <-P_Gneighbour %>% 
  inner_join(symphonia_stage, by = "ID")
 
indices <- indices %>% rbind(G_neighbour %>% dplyr::select(n_parcelle, n_carre, n_arbre, GNeighbors))
}



# Parcelle 16  ------------------------------------------------------------

# fusion coordinate subplot plot 16 ---------------------------------------

subplot <- t(matrix(c(1:25),5,5))

coordinate_raw  <- paracoutree %>% 
  filter(code_vivant ==1) %>% 
  select(n_parcelle, n_carre, n_arbre, X, Y, circ_corr, species) %>% 
  filter(n_parcelle == 16) %>%
  mutate(testcoord_x = 1, testcoord_y = 1)

for(i in 1:25){
  for( j in 1:nrow(subplot)){
    if(i %in% subplot[j,]){
      row <- j
      break
    }
  }
  for(k in 1:ncol(subplot)){
    if(i == subplot[row,k]){
      col = k
      break()
      
    }
  }
  increment_y = (5-row)*100
  increment_x = (col-1)*100
  
  #print(paste("subplot",i,"coordinates incremented by x:",increment_x,"y:",increment_y))
  
  coordinate_raw[which(coordinate_raw$n_carre == i), "testcoord_x"] <- coordinate_raw[which(coordinate_raw$n_carre == i), "X"] + increment_x 
  coordinate_raw[which(coordinate_raw$n_carre == i), "testcoord_y"] <- coordinate_raw[which(coordinate_raw$n_carre == i), "Y"] + increment_y
}



Paracoutest2 <- coordinate_raw %>%
  dplyr::select(n_parcelle, n_carre, n_arbre, testcoord_x, testcoord_y, circ_corr, species) %>% 
  rename(X = "testcoord_x",Y = "testcoord_y")


indices2 = matrix(ncol = 4)
indices2 <- as.data.frame(indices2); names(indices2) <- c("n_parcelle","n_carre", "n_arbre","GNeighbors")


  
Paracoutest2 <- Paracoutest2 %>% 
  mutate( PointWeight = (circ_corr)^2/(4*pi)) %>% 
  mutate( PointName = paste(n_parcelle, n_carre, n_arbre, sep="_")) %>%
  rename( PointType = species) %>% 
  na.omit()
rownames(Paracoutest2) <- Paracoutest2$PointName

#duplicated point verification
Paracoutest2 %>% select(X,Y) %>% duplicated() %>% which() %>% length()
# 171 nb duplicated tree
Paracoutest2 %>% select(X,Y, PointName) %>% duplicated() %>% which() %>% length()
# 0 it different tree probably really close

Paracoutest2 <- dbmss::wmppp(Paracoutest2, window = owin(c(min(Paracoutest2$X),max(Paracoutest2$X)), c(min(Paracoutest2$Y),max(Paracoutest2$Y)),
                                                       unitname=c("metres", "metres")))

# Définition du voisinage

rVoisinage <- 25


# 1 : sans correction des effets de bord
###########
# Matrice de distance
Distances <- pairdist(Paracoutest2)
# Matrice de voisinage (25m)
Voisins <- (Distances <= rVoisinage) # inférieur ou égale


# Elimination du point lui-même
diag(Voisins) <- FALSE

# Surface terrière des voisins, par colonne

Gvoisins <- apply(Voisins, 2, function(EstVoisin) sum(Paracoutest2[EstVoisin]$marks$PointWeight))

# 2 : correction des effets de bord
###########

# Facteur de correction
Correction <- function(NumPoint) {
  # Disque de 25m de rayon autour du point
  disc(radius=rVoisinage, centre=c(Paracoutest2$x[NumPoint], Paracoutest2$y[NumPoint])) %>% 
    # Intersection avec la parcelle
    intersect.owin(Paracoutest2$window) %>% 
    # Calul de la surface
    area -> VoisinageDansParcelle
  # Retour du facteur de correction
  return(pi * rVoisinage^2 / VoisinageDansParcelle)
}

Corrections <- vapply(1:Paracoutest2$n, Correction, 0)

# Surface terrière corrigée ajoutée au wmppp

Paracoutest2$marks$GNeighbors <- Gvoisins*Corrections

# Nombre de voisins
#Paracoutest$marks$nNeighbors <- colSums(Voisins)*Corrections


# # Carte
# plot(density(Paracoutest, weights = Paracoutest$marks$GNeighbors), main="Surface terrière des voisins")
# plot(Paracoutest[grep("Symphonia_", Paracoutest$marks$PointType)], which.marks = "PointWeight", add=TRUE)
# 
# # Surface terrière
# plot(density(Paracoutest, weights = Paracoutest$marks$PointWeight), main="Surface terrière")

# extraction des données 

P_Gneighbour2 <- as.data.frame(Paracoutest2) 
P_Gneighbour2 <- P_Gneighbour2 %>% 
  mutate(ID = rownames(P_Gneighbour2)) 
G_neighbour2 <-P_Gneighbour2 %>% 
  inner_join(symphonia_stage, by = "ID") # in environment only the last one with no symphonia

indices2 <- indices2 %>% rbind(G_neighbour2 %>% select(n_parcelle, n_carre, n_arbre, GNeighbors))


lBAL<- rbind(indices, indices2) %>% filter(!is.na(n_parcelle))

write.csv(lBAL, file ="C:\\Users\\emduc\\Desktop\\Drive\\symphostage\\data_merging\\Full_lBAL.csv" )




rm("G_neighbour", "P_Gneighbour","G_neighbour2", "P_Gneighbour2", "Distances", "Voisins", "i","j", "rVoisinage", "Gvoisins", "Corrections", "Paracoutest", "Paracoutest2")

# Hegyi -------------------------------------------------------------------



indices3 = matrix(ncol = 4)
indices3 <- as.data.frame(indices3); names(indices3) <- c("n_parcelle","n_carre", "n_arbre","GNeighbors")

for (i in 1:15) {
  
  Paracoutest <- paracoutree %>% 
    filter(code_vivant ==1) %>% 
    select(n_parcelle, n_carre, n_arbre, X, Y, circ_corr, species) %>% 
    filter(n_parcelle == i) %>% 
    mutate( PointWeight = (circ_corr)^2/(4*pi)) %>% 
    mutate( PointName = paste(n_parcelle, n_carre, n_arbre, sep="_")) %>%
    rename( PointType = species) %>% 
    na.omit()
  rownames(Paracoutest) <- Paracoutest$PointName
  Paracoutest <- dbmss::wmppp(Paracoutest, window = owin(c(0,250), c(0,250),
                                                         unitname=c("metres", "metres")))
  
  
  
  # Définition du voisinage
  
  rVoisinage <- 25
  
  
  # 1 : sans correction des effets de bord
  ###########
  # Matrice de distance
  Distances <- pairdist(Paracoutest)
  # Matrice de voisinage (25m)
  Voisins <- (Distances <= rVoisinage) # inférieur ou égale
  
  
  # # Variante : seulement les gros arbres
  # Voisins <- (Distances <= rVoisinage) & (Paracoutest$marks$PointWeight > 1000)
  # 
  # # Variante 2 : les symphonia
  # EstSymphonia <- logical(Paracoutest$n) #vecteur de la longueur de paracoutest ici 3540
  # EstSymphonia[grep("Symphonia_", Paracoutest$marks$PointType)] <- TRUE #Èremplacement par TRUE des place de dans le vecteur qui
  # #corresponde à celle dans lesquelles il y a le pattern Symphonia_
  # Voisins <- (Distances <= rVoisinage) & EstSymphonia # on selectionne dans la matrice voisin ceux qui répondent aux 2 conditions
  
  
  # Elimination du point lui-même
  diag(Voisins) <- FALSE
  
  # Surface terrière des voisins, par colonne
  
  Hegyi <- sapply(1:ncol(Voisins), function(numPoint) sum(Paracoutest[Voisins[, numPoint]]$marks$PointWeight/Distances[Voisins[, numPoint], numPoint]))
  
  # 2 : correction des effets de bord
  ###########
  
  # Facteur de correction
  Correction <- function(NumPoint) {
    # Disque de 25m de rayon autour du point
    disc(radius=rVoisinage, centre=c(Paracoutest$x[NumPoint], Paracoutest$y[NumPoint])) %>% 
      # Intersection avec la parcelle
      intersect.owin(Paracoutest$window) %>% 
      # Calul de la surface
      area -> VoisinageDansParcelle
    # Retour du facteur de correction
    return(pi * rVoisinage^2 / VoisinageDansParcelle)
  }
  
  Corrections <- vapply(1:Paracoutest$n, Correction, 0)
  
  # Surface terrière corrigée ajoutée au wmppp

  
Paracoutest$marks$GNeighbors <- Hegyi*Corrections
  
  # Nombre de voisins
  #Paracoutest$marks$nNeighbors <- colSums(Voisins)*Corrections
  
  
  # # Carte
  # plot(density(Paracoutest, weights = Paracoutest$marks$GNeighbors), main="Surface terrière des voisins")
  # plot(Paracoutest[grep("Symphonia_", Paracoutest$marks$PointType)], which.marks = "PointWeight", add=TRUE)
  # 
  # # Surface terrière
  # plot(density(Paracoutest, weights = Paracoutest$marks$PointWeight), main="Surface terrière")
  
  # extraction des données 
  
  P_Gneighbour <- as.data.frame(Paracoutest) 
  P_Gneighbour <- P_Gneighbour %>% 
    mutate(ID = rownames(P_Gneighbour)) 
  G_neighbour <-P_Gneighbour %>% 
    inner_join(symphonia_stage, by = "ID")
  
  indices3 <- indices3 %>% rbind(G_neighbour %>% dplyr::select(n_parcelle, n_carre, n_arbre, GNeighbors))
}




# Parcelle 16 -------------------------------------------------------------



indices4 = matrix(ncol = 4)
indices4 <- as.data.frame(indices4); names(indices4) <- c("n_parcelle","n_carre", "n_arbre","GNeighbors")


for (j in 1:25) {
  
  Paracoutest2 <- paracoutree %>% 
    filter(code_vivant ==1) %>% 
    select(n_parcelle, n_carre, n_arbre, X, Y, circ_corr, species) %>% 
    filter(n_parcelle == 16) %>% 
    filter(n_carre == j) %>% 
    mutate( PointWeight = (circ_corr)^2/(4*pi)) %>% 
    mutate( PointName = paste(n_parcelle, n_carre, n_arbre, sep="_")) %>%
    rename( PointType = species) %>% 
    na.omit()
  rownames(Paracoutest2) <- Paracoutest2$PointName
  Paracoutest2 <- dbmss::wmppp(Paracoutest2, window = owin(c(min(Paracoutest2$X),max(Paracoutest2$X)), c(min(Paracoutest2$Y),max(Paracoutest2$Y)),
                                                         unitname=c("metres", "metres")))
  
  
  
  
  # Définition du voisinage
  
  rVoisinage <- 25
  
  
  # 1 : sans correction des effets de bord
  ###########
  # Matrice de distance
  Distances <- pairdist(Paracoutest2)
  # Matrice de voisinage (25m)
  Voisins <- (Distances <= rVoisinage) # inférieur ou égale
  
  
  # # Variante : seulement les gros arbres
  # Voisins <- (Distances <= rVoisinage) & (Paracoutest2$marks$PointWeight > 1000)
  # 
  # # Variante 2 : les symphonia
  # EstSymphonia <- logical(Paracoutest2$n) #vecteur de la longueur de Paracoutest2 ici 3540
  # EstSymphonia[grep("Symphonia_", Paracoutest2$marks$PointType)] <- TRUE #Èremplacement par TRUE des place de dans le vecteur qui
  # #corresponde à celle dans lesquelles il y a le pattern Symphonia_
  # Voisins <- (Distances <= rVoisinage) & EstSymphonia # on selectionne dans la matrice voisin ceux qui répondent aux 2 conditions
  
  
  # Elimination du point lui-même
  diag(Voisins) <- FALSE
  
  # Surface terrière des voisins, par colonne
  

Hegyi <- sapply(1:ncol(Voisins), function(numPoint) sum(Paracoutest2[Voisins[, numPoint]]$marks$PointWeight/Distances[Voisins[, numPoint], numPoint]))
  
  # 2 : correction des effets de bord
  ###########
  
  # Facteur de correction
  Correction <- function(NumPoint) {
    # Disque de 25m de rayon autour du point
    disc(radius=rVoisinage, centre=c(Paracoutest2$x[NumPoint], Paracoutest2$y[NumPoint])) %>% 
      # Intersection avec la parcelle
      intersect.owin(Paracoutest2$window) %>% 
      # Calul de la surface
      area -> VoisinageDansParcelle
    # Retour du facteur de correction
    return(pi * rVoisinage^2 / VoisinageDansParcelle)
  }
  
  Corrections <- vapply(1:Paracoutest2$n, Correction, 0)
  
  # Surface terrière corrigée ajoutée au wmppp
  
  
  Paracoutest2$marks$GNeighbors <- Hegyi*Corrections
  
  # Nombre de voisins
  #Paracoutest2$marks$nNeighbors <- colSums(Voisins)*Corrections
  
  
  # # Carte
  # plot(density(Paracoutest2, weights = Paracoutest2$marks$GNeighbors), main="Surface terrière des voisins")
  # plot(Paracoutest2[grep("Symphonia_", Paracoutest2$marks$PointType)], which.marks = "PointWeight", add=TRUE)
  # 
  # # Surface terrière
  # plot(density(Paracoutest2, weights = Paracoutest2$marks$PointWeight), main="Surface terrière")
  
  # extraction des données 
  
  P_Gneighbour2 <- as.data.frame(Paracoutest2) 
  P_Gneighbour2 <- P_Gneighbour2 %>% 
    mutate(ID = rownames(P_Gneighbour2)) 
  G_neighbour2 <-P_Gneighbour2 %>% 
    inner_join(symphonia_stage, by = "ID")
  
  indices4 <- indices4 %>% rbind(G_neighbour2 %>% select(n_parcelle, n_carre, n_arbre, GNeighbors))
}

Competition <- rbind(indices3, indices4) %>% filter(!is.na(n_parcelle))

write.csv(Competition, file ="C:\\Users\\emduc\\Desktop\\Drive\\symphostage\\data_merging\\Full_Hegyi.csv" )


