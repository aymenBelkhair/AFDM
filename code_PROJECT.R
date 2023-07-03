heart <- read.csv("heart_prediction.csv",sep = ",")

summary(heart)
getwd()



############## PHASE 1 : Prétraitement des données ##############

###### ACP sur les variables quantitatives

#fonction pour centrage et réduction 

CR <- function(x) {
  n <- length(x)
  m <- mean(x)
  v <- (n-1)/n*var(x)
  return((x-m)/sqrt(v))
}


#appliquer la fonction sur les variables continues

Vquanti <- data.frame(heart$ï..age,heart$trestbps,heart$chol,heart$thalach,heart$oldpeak)
print(head(Vquanti))
heart.cont <- data.frame(lapply(Vquanti,CR))
print(heart.cont)



###### ACM sur les variables qualitatives


#codage disjonctif complet

Vquali <- data.frame(heart$sex,heart$cp,heart$fbs,heart$restecg,heart$exang,
                     heart$slope,heart$ca,heart$thal)

library(ade4)

heart.disc <- acm.disjonctif(Vquali)

#fonction pour pondération des indicatrices
PF <- function (x){
  m <- mean(x)
  return(x/sqrt(m))
}

#appliquer la pondération sur les indicatrices
heart.disc.pond <- data.frame(lapply(heart.disc,PF))

print(heart.disc.pond)






#données transformées envoyées à l'AFD

heart_sans_target <- cbind(heart.cont,heart.disc.pond)

heart.pour.afd <- cbind(heart.cont,heart.disc.pond,heart$target)

df <- data.frame(round(heart.pour.afd,3))

dp <- data.frame(heart.disc.pond)



rownames(heart.pour.afd) <- row.names(heart)
print(round(heart.pour.afd,3))


#### Matrice pré-traité pour AFD

df <- data.frame(round(heart.pour.afd,3))


write.csv(df, file = "C:/Users/user/Documents/file3.csv", row.names = T)

###### les données centrées pour enveoyé a l'ACP

XX= scale(df[,-31],scale = F)

dfx <- data.frame(XX)

# Chargement du package
library(MASS)

# Effectuer l'analyse factorielle discriminante
res.afd3=lda(df$heart.target~.,dfx)


# Afficher les résultats
print(res.afd3)
plot(res.afd3)



################### Comparaison ##############


# Création du modèle d'analyse factorielle discriminante

modele1 <- lda(heart$target ~ ., data = heart)


# Affichage des résultats

print(modele1)

plot(modele1)




Vquanti <- data.frame(heart$age,heart$trestbps,heart$chol,heart$thalach,heart$oldpeak)
Vquali <- data.frame(heart$sex,heart$cp,heart$fbs,heart$restecg,heart$exang,heart$slope,heart$ca,heart$thal)




# calcule de la matrice G 

G11=mean(heart.pour.afd$heart.ï..age[heart.pour.afd$`heart$target`=="0"])
G21=mean(heart.pour.afd$heart.ï..age[heart.pour.afd$`heart$target`=="1"])
print(c(G11,G21))

# les memes valeurs trouver dans package donc on extrait la matrice G A partir du résultat du lda 

group_means <- res.afd3$means
G <- group_means
# Print the group means matrix
print(G)




x <- heart_sans_target



# la matrice de variances inter-classes B = tG*??*G d'ou  :

n1 <- sum(heart.pour.afd$`heart$target`==1)
n2 <- sum(heart.pour.afd$`heart$target`==0)
n<- n1+n2

# delta

delta <- matrix(c(n1/n, 0, 0, n2/n), nrow = 2, ncol = 2)

B = (t(G)%*%delta%*%G)


#la matrice de variance total V:
TX <- t(XX)
V1=(1/n)*TX%*%XX

V1  

# inverse de la matrice V

VI = ginv(V1)
M = VI%*%B


# les valeurss propores :
val_p = Re(eigen(M)$values)  #Re pour eliminer la partie imaginaire
vp1=Re(eigen(M)$values[1])
# vp2=Re(eigen(M)$values[2]) puisque la 2eme valeur propre contient une valeur
# négligable on travaille seulement avec la 1er valeur propre 



# Part d'inertie des valeurs propres

Sum_val = sum(val_p)
part_inertie = (vp1)/Sum_val  #100 





# les vecteurs propores :

Vect_p =Re(eigen(M)$vectors)

vec_p_1 =Re(eigen(M)$vectors[,1])
#vec_p_2 =Re(eigen(M)$vectors[,2])
vec_p_1


#les vecteurs principaux : 

norme_vect_propre_1= sqrt(t(vec_p_1)%*%VI%*%vec_p_1)
#norme_vect_propre_2= sqrt(t(vec_p_2)%*%VI%*%vec_p_2) 


vectprin1 = (1/norme_vect_propre_1)*vec_p_1
#vectprin2 = (1/norme_vect_propre_2)*vec_p_2

vectprin1
#les coordonées des centres de gravités (Composantes Principales):

C11=G%*%VI%*%vectprin1  # ON travaille seulement avec C11 puisque lambda1 = 100% D'inertie initiale
#C22=G%*%VI%*%vectprin2

C11



plot(C11)

#les coordonées des individus :

CI1=XX%*%VI%*%vectprin1
#CI2=XX%*%VI%*%vectprin2

print(CI1)
#====================
library(ggplot2)

# Création d'un dataframe avec les coordonnées C11[1], C11[0] et CI1
df <- data.frame(x = c(C11[1], C11[2], CI1), groupe = c(1, 0, rep("Autres", length(CI1))))

# Création du graphique de dispersion avec une couleur différente pour chaque groupe
ggplot(df, aes(x = x, y = 0, color = groupe)) +
  geom_point(size = ifelse(df$groupe %in% c(1, 0), 5, 2)) +
  labs(x = "Coordonnées", y = "") +
  ggtitle("Représentation graphique des coordonnées C11[1], C11[0] et CI1") +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme(legend.position = "none")



plot(CI1, pch = 16, col = df$heart.target)  # Remplacez 'y' par votre vecteur de classes
text(CI1, labels = df$heart.target, pos = 2) 
new_points <- c(C11[1], C11[2])   # Coordonnées x des nouveaux points

points(new_points,  pch = 16, col = "blue")  # Ajouter les points bleus
text(new_points, label,  pos = 2)  # Ajouter les étiquettes


plot(CI1, pch = 16, col = df$heart.target)
text(CI1, labels = df$heart.target, pos = 2)
new_points <- c(C11[1], C11[2])   # Coordonnées x des nouveaux points

points(new_points, pch = c(16, 16), cex = c(2, 2), col = c("red", "green"), font = 2)  # Ajouter les points en rouge et vert de taille plus grande et en gras
text(new_points, label, pos = 2)
#----------------------------------------------

library(ggplot2)

points <- data.frame(x = CI1, labels = "")

ggplot(points, aes(x = x, y = 0, z = 0)) +
  geom_point(size = 2, aes(color = df$heart.target), shape = 16) +
  geom_text(aes(label = labels), vjust = -1) +
  labs(x = "Valeur de CI1", y = "", z = "") +
  ggtitle("Représentation des individus sur un axe CI1") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()
#------------------------------------------
plot(CI1, rep(0, length(CI1)), pch = 16, col = df$heart.target, xlim = c(-4, 4), ylim = c(0, 1), ylab = "Corrd individus", xlab = "Corrd centre de graviter")
text(CI1, rep(0, length(CI1)), labels = df$heart.target, pos = 2)
new_points <- c(C11[1], C11[2])   # Coordonnées x des nouveaux points

points(new_points, rep(0, length(new_points)), pch = c(16, 16), cex = c(2, 2), col = c("red", "green"), font = 2)
text(new_points, 0, label, pos = 2)


###### exemple sur un seul axe


library(ggplot2)

points <- data.frame(x=c(CI1),labels = "")

ggplot(points, aes(x = x, y = 0,z=0)) +
  geom_point() +
  geom_text(aes(label = labels), vjust = -1) +
  labs(x = "Valeur de x", y = "",z="") +
  ggtitle("Ensemble de points sur un axe C1")


########## Affectation d'un nouveau individu


# l'objectif de cette phase est de prédire ......


########## Affectation d'un nouveau individu






x=heart_sans_target


########## Affectation d'un nouveau individu


# l'objectif de cette phase est de prédire l'appartenance d'un nouveau individu , on connait leur variable explicatif ms pas leur classe



# prétraitement de l'observavtion i

s <- c(45,0,2,120,300,1,1,140,1,3,1,0,2)
header <- c("age","sex","cp","trestbps","chol","fbs","restecg","thalach",
            "exang","oldpeak","slope","ca","thal")
names(s) <- header
dss <- data.frame(t(s))
dss



DSS_quanti <- data.frame(dss$age,dss$trestbps,dss$chol,dss$thalach,dss$oldpeak)
DSS_quali <- data.frame(dss$sex,dss$cp,dss$fbs,dss$restecg,dss$exang,
                        dss$slope,dss$ca,dss$thal)
dss1 <- cbind(DSS_quanti,DSS_quali)
dss1_center <-c(45-mean(x$heart.ï..age),120-mean(x$heart.trestbps),
                300-mean(x$heart.chol),140-mean(x$heart.thalach),
                3-mean(x$heart.oldpeak),
                1-mean(x$heart.sex.0) ,0 ,    #sex
                0,0-mean(x$heart.cp.1),1-mean(x$heart.cp.2),0,   #cp
                0,1-mean(x$heart.fbs.1),        #fbs
                0,1-mean(x$heart.restecg.1),0,       #restecg
                0,1-mean(x$heart.exang.1),          #exang
                0,1-mean(x$heart.slope.1),0,         #slope
                1-mean(x$heart.ca.0),0,0,0,0,0,0,1-mean(x$heart.thal.2),0 )#thal ca
axe11=dss1_center%*%VI%*%vectprin1


#### en deduit que l'individu i est attiré par le groupe 0 (L'absence de maladie cardiaque)


###########