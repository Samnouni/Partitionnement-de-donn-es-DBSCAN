#-----------------------------------------------
# chargement des packages et fonctions utilises 
#----------------------------------------------

library(nnet)
library(ggplot2)
library("FactoMineR")
library("factoextra")
require(jpeg)

#----------------------------------------------------------------------------------------------------------------
                                          # Exercice 1 
#----------------------------------------------------------------------------------------------------------------

#---------------------------------
# chargement du fichier de donnees
#--------------------------------- 

animal <- read.csv("/home/samnouni/Bureau/M2 ISN/methode daprentissage/projet/animaux_I.csv",sep=",", header=TRUE) 

#---------------------------------
# Variable en factor
#--------------------------------- 

for (i in 2:ncol(animal)){animal[,i]=as.factor(animal[,i])} 

#------------------------------------
# Recodage de la variable class_type
#------------------------------------

levels(animal$class_type)=c( "mammifère","oiseau","reptile"," poisson","amphibien","insecte","invertébré ")
animal=animal[,-1]

#------------------------------------
# refaire une analyse descriptive
#------------------------------------

str(animal)
summary(animal)
head(animal)

#------------------------------------
# Calcul errreur 
#------------------------------------

tx_er=function(pred,vrais)
{
  mc=table(pred,vrais) 
  1-sum(diag(mc))/sum(mc)
} 

#-----------------------------------------
# Calcule Taux d'erreur du modèle initial
#-----------------------------------------

Erreurs=function(animal){
  N <- nrow(animal)
  ntrain <- floor(0.80*N) # partie entier 
  ntest <- N - ntrain

  indices <- sample(1:N,ntrain,replace = FALSE)

  animal_train <- animal[indices,]
  animal_test <- animal[-indices,]
  log.animal= multinom(class_type~.,data=animal_train)
  predrf <- predict(log.animal,newdata=animal_test,type="class")
  te_rf <- tx_er(predrf,animal_test$class_type)
  return(te_rf)
}

#---------------------------------------------
# Graphique pour l'erreur en chaque itération
#---------------------------------------------

Erreurs_modele_initial=sapply(1:50 ,FUN=function(i){Erreurs(animal)})
plot(100*Erreurs_modele_initial,main="Le taux d'erreurs du modèle initial",ylab="Taux d'erreur %",xlab="Itération")
mean(Erreurs_modele_initial)

#---------------------------------------------
# Choix d'un modèle avec moins de variables
#---------------------------------------------

log.animal= multinom(class_type~.,data=animal_train)
summary(log.animal) # Pour afficher les resultats du modele
stepAIC(log.animal, trace = TRUE, data = animal_train) # pour choisir le modéle le plus pertinent au snes du critère d'AIC

#---------------------------------------------
# Premier modèle trouvé
#---------------------------------------------

log.animal1= multinom(class_type~feathers+milk+backbone+breathes,data=animal_train)
predrf1 <- predict(log.animal1,newdata=animal_test,type="class")
te_rf1 <- tx_er(predrf1,animal_test$class_type)
te_rf1

#---------------------------------------------
# Deuxième modèle trouvé 
#---------------------------------------------

log.animal2= multinom(class_type~hair+feathers+aquatic+breathes+venomous+tail,data=animal_train)
predrf2 <- predict(log.animal2,newdata=animal_test,type="class")
te_rf2 <- tx_er(predrf2,animal_test$class_type)
te_rf2

#---------------------------------------------
# ACM (Étude de la liaison entre les variables)
#---------------------------------------------

animal=animal[,-18]
res.mca <- MCA (animal, graph = FALSE)
eig.val <- get_eigenvalue(res.mca) # Valeur propre
fviz_screeplot (res.mca, addlabels = TRUE, ylim = c (0, 45))
fviz_mca_var (res.mca, choice = "mca.cor",repel = TRUE,ggtheme = theme_minimal (),axes = c(3,2))
fviz_mca_var (res.mca,repel = TRUE,ggtheme = theme_minimal (),head=T,axes = c(1,2))

#---------------------------------------------
# Modèle final
#---------------------------------------------

log.animal_final1= multinom(class_type~eggs+aquatic+legs+toothed+catsize+predator,data=animal_train)
predrf_final1 <- predict(log.animal_final1,newdata=animal_test,type="class")
te_rf_final1 <- tx_er(predrf_final1,animal_test$class_type)
te_rf_final1
stepAIC(log.animal_final1, trace = TRUE, data = animal_train) # pour choisir le modéle le plus pertinent au snes du critère d'AIC

#---------------------------------------------
#Fonction pour calculer l'erreur
#---------------------------------------------

Erreurs_modele_final=function(animal){
  N <- nrow(animal)
  ntrain <- floor(0.80*N) # partie entier 
  ntest <- N - ntrain
  
  indices <- sample(1:N,ntrain,replace = FALSE)
  
  animal_train <- animal[indices,]
  animal_test <- animal[-indices,]
  log.animal= multinom(class_type~eggs+aquatic+legs+toothed+catsize,data=animal_train)
  predrf <- predict(log.animal,newdata=animal_test,type="class")
  te_rf <- tx_er(predrf,animal_test$class_type)
  return(te_rf)
}

Erreurs_modele_final=sapply(1:50 ,FUN=function(i){Erreurs_modele_final(animal)})
plot(100*Erreurs_modele_final,main="Le taux d'erreurs du modèle final",ylab="Taux d'erreur %",xlab="Itération")
mean(Erreurs_modele_final)


#### Partie 2 ####

data = data.frame() 
# On récupère les noms des fichiers
setwd("~/Bureau/Méthodes_Dapprentissag/Projet")
animaux = list.files(path="animaux_II",full.names = F, recursive = F)
# choisir en hasard 2004 images
indices_image =sample(1:length(animaux),0.423*length(animaux))
names_animaux=animaux[indices_image]
nbimages=length(names_animaux)

for (i in 1:nbimages)
{
  files=paste0("animaux_II/",names_animaux[i])
  img=readJPEG(files)
  # transformation de chaque image en vecteur
  var1=as.vector(img[,,1])
  var2=as.vector(img[,,2])
  var3=as.vector(img[,,3])
  v=c(var1,var2,var3)
  v=as.data.frame(t(v))
  
  data =rbind(data,v )# on combine les vecteurs
  
  
}
# on renomme chaque ligne par nom de fichier
rownames(data)=names_animaux


# Application de l'acp sur notre base
acp = prcomp(data, scale. = F, center = F, retx = TRUE)

# on garde les axes dont les valeurs propres sont supérieures à la moyenne
nbaxe = sum(acp$sdev**2 >= mean(acp$sdev**2))
pourcentage_inertie =sum(acp$sdev[1:nbaxe]**2)/sum(acp$sdev**2) # Pourcentage d'inertie
pourcentage_inertie


source("~/Bureau/Méthodes_Dapprentissag/DBSCAN/color_utils.R")
# On teste la performance de l'Acp en fonction de nombre d'axes
img_reelle = array(unlist(data[1,]),dim=c(128,128,3))
display(img_reelle)

# les coordonnées des images dans la base formée par les 44 premier axes
data_acp = as.data.frame(acp$x[,1:nbaxe])
# Transformation à la base initiale
mat_pass = acp$rotation[,1:nbaxe] # Matrice de passage
images_reduite <- as.matrix(data_acp) %*% t(mat_pass) 
img_reduite <- array(unlist(images_reduite[1,]),dim=c(128,128,3))
img_reduite[img_reduite<0] <- 0
img_reduite[img_reduite>1] <- 1
display(img_reduite)

# les coordonnées des images dans la base formée par les 200 premiers axes
nbaxe=200
data_acp = as.data.frame(acp$x[,1:nbaxe])
# Transformation à la base initiale
mat_pass = acp$rotation[,1:nbaxe] # Matrice de passage
images_reduite <- as.matrix(data_acp) %*% t(mat_pass) 
img_reduite <- array(unlist(images_reduite[1,]),dim=c(128,128,3))
img_reduite[img_reduite<0] <- 0
img_reduite[img_reduite>1] <- 1
display(img_reduite)

### On va travailler sur les 200 premiers axes ####
data_acp=as.data.frame(acp$x[,1:nbaxe])
data_acp_plan1=data_acp[,1:2] # La projection de données sur le premier plan
colnames(data_acp_plan1)=c("x1","x2")
### DBSCAN ###

dKnn=kNNdist(data_acp,k=10)

plot(sort(dKnn,decreasing = T),type = 'l',xlim = c(0,500)) 
dev.off()
eps=54
dbing=dbscan(data_acp,eps=eps,minPts = 10)
clust=dbing$cluster
clust[clust==0]="bruit"
data_acp_plan1_DBSCAN=cbind(data_acp_plan1,clust)
ggplot(data = data_acp_plan1_DBSCAN,aes(x1,x2,color=as.factor(clust))) + geom_point()+labs(color = "Classe")


### CAH ###
dist_d1 <- dist(data_acp,method = "euclidean")
hclust_1 <- hclust(dist_d1,method="ward.D")
plot(hclust_1)
# on garde 5 classes
clusters <- rect.hclust(hclust_1, k=5)
cah_clust=rep(1,nrow(data_acp_plan1))
cah_clust[clusters[[2]]]=2
cah_clust[clusters[[3]]]=3
cah_clust[clusters[[4]]]=4
cah_clust[clusters[[5]]]=5
data_acp_plan1_CAH=cbind(data_acp_plan1,cah_clust)
ggplot(data = data_acp_plan1_CAH,aes(x1,x2,color=as.factor(cah_clust))) + geom_point()+labs(color = "Classe")

### K-NN ###
sse <- numeric()
for(i in 1:14){
  sse[i]=sum(kmeans(data_acp,centers=i,nstart = 2)$withinss)
}
plot(1:14,sse,pch=20,type="b",xlab = "Nombres de classes",ylab = "Inertie intra-classes")
# on choisit par exemple 4 ou 5
clust_KNN =kmeans(data_acp,centers=5, nstart = 2)
data_acp_plan1_KNN = cbind(data_acp_plan1,clust=clust_KNN$cluster)
ggplot(data = data_acp_plan1_KNN,aes(x1,x2,color=as.factor(clust))) + geom_point()+labs(color = "Classe")


