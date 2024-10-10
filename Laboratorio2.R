#laboratorio2.R
library(tidyverse)
data(starwars)

starwars %>% select(-name)
starwars %>% select(contains("_"))
starwars %>% select(starts_with("s"))
homeworld <- starwars %>% select(name, homeworld)

human <- starwars %>% filter(species == "Human")
starwars %>% filter(species == "Human", homeworld == "Tatooine")
starwars_nandroids <- starwars %>% filter(species != "Droid")
#77 registros cumplen las condiciones finales

starwars %>% group_by(species) %>% tally()
starwars %>% group_by(species, gender) %>% tally()
table_gender <- starwars %>% group_by(species, gender) %>% tally()

starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T), mean_mass = mean(mass, na.rm = T))
#desviacion estandar: 
starwars %>% group_by(species) %>% summarise(sd_height = sd(height, na.rm = T), sd_mass = sd(mass, na.rm = T))

ggplot(starwars, aes(height, mass)) + geom_point()
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()

#para el primer gráfico eliminando el personaje con masa muy grande 
starwars_good <- starwars %>% filter(name != "Jabba Desilijic Tiure")
ggplot(starwars_good, aes(height, mass)) + geom_point(colour = "pink")
#el gráfico esta en el pdf

toy <- read_csv("C:/Users/Minerva/Downloads/toy.csv")

#agrupación por sexo
toy %>% group_by(Sex) %>% summarise(mean_weight = mean(Weight_Kg, na.rm = T), mean_height = mean(Height_cm, na.rm = T), mean_IMC = mean(IMC, na.rm = T), mean_IAS = mean(IAS, na.rm = T), mean_CCintura = mean(Ccintura, na.rm = T))

tabla_fem <- toy %>% filter(Sex == "Women") %>% summarise(mean_weight = mean(Weight_Kg, na.rm = T), mean_height = mean(Height_cm, na.rm = T), mean_IMC = mean(IMC, na.rm = T), mean_IAS = mean(IAS, na.rm = T), mean_CCintura = mean(Ccintura, na.rm = T))
tabla_women <- toy %>% filter(Sex == "Women")
tabla_women %>% tally()
#58 registros cumplen las condiciones
tabla_women %>% filter(IMC_clas == "Overweight") %>% tally()
#9 de ellos tienen sobrepeso

ggplot(toy, aes(IMC, Weight_Kg)) + geom_point(colour="blue")
ggplot(toy %>% filter(IMC_clas == c("Overweight", "Obesity")), aes(IMC, Weight_Kg)) + geom_point(colour="blue")
#gráficos en el pdf

install.packages("ape")
install.packages("phangorn")
install.packages("phytools")