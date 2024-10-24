library(ape)
library(phangorn)
library(phytools)

fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", format = "FASTA", type = "AA")
fraxatin

matrizdist <- as.AAbin(fraxatin) #fraxatin en matriz como aminoacidos
matrizdist <- dist.aa(matrizdist) #matriz mide distancias
matrizdist

#los numeros que me da la matriz son similitudes en aas 

arbolUPGMA <- upgma(matrizdist) #bajo funcion upgma con info de matrizdist meto en arbolUPGMA
plot(arbolUPGMA)
#Si la longitud de dos ramas es indéntica, significa que 
#las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0

#metodo NJ (2) union de vecinos
arbolNJ <- nj(matrizdist)
plot(arbolNJ)

#diseño propio
plot(arbolUPGMA, type="p", cex=0.6, edge.width = 2, edge.color="orange", font=3)
#type: tipo de lineas, cex: letra tamano, edgew: grosor, edgeclor, font: tipo de letra

plotTree(arbolNJ) #phytool plot
plotTree(ladderize(arbolNJ)) #organizar por grupos mas comunes

write.tree(arbolNJ, file= "file_tree.nex") #guarda arbol
read.tree(file= "file_tree.nex") #leer lo que dice el arbol info

arbolNJraiz <- root(arbolNJ, outgroup= "Ornitorrinco", r= TRUE) #enraizar/comparar
plot(arbolNJraiz)

arbolUPGMAraiz <- root(arbolUPGMA, outgroup = "Ornitorrinco", r= TRUE)
plot(arbolUPGMAraiz)

layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1)) #par: parametros, mar: margenes 
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex = 0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex = 0.4)

#ARBOL PARSIMONIA:disminuir el número de pasos que explican un árbol evolutivo
#contando el número de cambios de cada uno de los caractere 
#El árbol con menor número de cambios es el ideal.
#no cuenta aas repertitivos o variables

parsimony(arbolUPGMAraiz, fraxatin)
parsimony(arbolUPGMA, fraxatin)
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)

fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE) #CANTIDAD CON MENOR LONGITUD
fraxatin_parsimonia #igual longitud

fraxatin_parsimoniaR <-root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex= 0.8)


estrictode100 <- consensus(fraxatin_parsimoniaR, p=1)
plot(estrictode100, cex=.6)

estrictode50 <- consensus(fraxatin_parsimoniaR, p=0.5)
plot(estrictode50, cex=.6) #intento menor probabilidad

estrictode10 <- consensus(fraxatin_parsimoniaR, p=0.1)
plot(estrictode10, cex=.6) #intento menor probabilidad

arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs=10)
#bs: numero de replicas, fun: funcion, pratchet: funcion de parsimony tree que haga menor numero

plot(arbolesbootstrap, cex=.7)
estricto60 <- consensus(arbolesbootstrap, p=0.6)
plot(estricto60, cex=.6)

arbolazar <- rtree(n=11, tip.label = names(fraxatin))
plot(arbolazar, cex=.5)
#rtree random tree

arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex=.5); add.scale.bar()

ajustado <- pml(arbolazarR, fraxatin)
ajustado 
# dio -4527.608 lo que es la verosimilitud en este arbol
#optim.pml de phangorn es similitud y evolucion aa

ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")
ajustadoconDay$tree #encontrar el arbol oculto #dayhoff es matriz de sustitucion

ajustadoconDayraiz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraiz), cex=.5); add.scale.bar()

ajustadoconBlo <- optim.pml (object = ajustado, model = "Blosum62", rearrangement = "ratchet")
ajustadoconJTT <- optim.pml(object= ajustado, model="JTT", rearrangement = "ratchet")

#comparar los modelos AIC, queremos num bajo JTT mejor para evaluar distancia entre sec de prot y verosimilitud
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)

mejorarbol <- optim.pml(
  object = ajustadoconDay,
  model = "JTT",
  rearrangement = "ratchet") #no es la unica opcion de optim.pml 
mejorarbol

mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex=0.5); add.scale.bar()
#mejor arbol de maxima verosimilitud 
