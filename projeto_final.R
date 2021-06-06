#requisitando pacotes
require(raster)
require(rgdal)
#require(sp)
require(rgeos)
require(ggplot2)
require(dplyr)
require(reshape2)
require(randomForest)
require(e1071)
require(dplyr)
require(caret)
require(caTools)
require(prettymapr)


#carregar imagem das bandas
arquivos = list.files(path = "imagens_Sentinel_2A/imagens_mescladas/", pattern = "_B", full.names = T)
img = stack(arquivos)

#carregar shapefile
brasilia = readOGR("shapefiles/53UF2500G.shp")
plot(brasilia)

#verificar o sistema de coordenadas das imagens e do shapefile
crs(img)
crs(brasilia)

#alterar o sistema de coordenadas do shapefile
brasilia.utm = spTransform(x = brasilia, CRSobj = crs(img))
crs(brasilia.utm)


#recortar para area delimitada pelo shapefile
brasilia.mask = mask(x = img, mask = brasilia.utm)
brasilia.crop = crop(brasilia.mask, brasilia.utm)

#plotar as as imagens recortadas
par(mfrow = c(1,2))
plot(brasilia.mask$E_B8A)
plot(brasilia.crop$E_B8A)

#salvar a imagem recortada
writeRaster(x = brasilia.crop, filename = "arquivos_gerados/imagem_recortada.tif", overwrite = TRUE)

#_______________________________________________________________________________

#carregando imagem recortada salva
img = stack("arquivos_gerados/imagem_recortada.tif")

#definindo o nome das bandas
names(img) = c("B02", "B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12")
names(img)

#plotando a Banda 8A
plot(img$B8A)

#calculando os indices de vegetacao
img$NDVI = (img$B8A - img$B04) / (img$B8A + img$B04)
img$simple_ratio = img$B8A / img$B04
img$EVI = 2.5 *((img$B8A - img$B04)/10000) / (img$B8A / 10000 + 6 * img$B04/10000 - 7.5 * img$B02/10000 + 1)
img$NDWI = (img$B03 - img$B12) / (img$B03 + img$B12)

#plotando os indices
par(mfrow = c(2, 2))
plot(img$NDVI, col = gray(0:100 / 100), main = "NDVI")
plot(img$simple_ratio, col = gray(0:100 / 100), main = "SR")
plot(img$EVI, col = gray(0:100 / 100), main = "EVI")
plot(img$NDWI, col = gray(0:100 / 100), main = "NDWI")

par(mfrow = c(1, 2))
plot(img$NDVI, col = gray(0:100 / 100), main = "NDVI")
plot(img$NDWI, col = gray(0:100 / 100), main = "NDWI")

#salvando imagem recortada com indices
writeRaster(x = img, filename = "arquivos_gerados/brasilia_com_indices.tif")

#salvando o nome das bandas
nomes_img = names(img)
write.csv(x = nomes_img, file = "arquivos_gerados/nomes_bandas.csv")

#_______________________________________________________________________________

#carregar imagem recortada com indices
brasilia_com_indices = stack("arquivos_gerados/brasilia_com_indices.tif")
names(brasilia_com_indices)

#alterar nome das bandas
nomes_bandas = read.table("arquivos_gerados/nomes_bandas.csv", header = T, sep = ',')
print(nomes_bandas)
names(brasilia_com_indices) = nomes_bandas[, 2]
names(brasilia_com_indices)

#importar arquivo shapefile com amostras
amostras = readOGR("amostras/amostras.shp")
View(data.frame(amostras))
unidos.shp = gUnaryUnion(spgeom = amostras, id = amostras$Classe)
unidos.shp
atributos = extract(x = brasilia_com_indices, y = unidos.shp)
names(unidos.shp)

#definir os atributos das amostras
Agricultura = data.frame(Classe = "Agricultura", atributos[[1]])
Agua = data.frame(Classe = "Água", atributos[[2]])
Area_Urbana = data.frame(Classe = "Área Urbana", atributos[[3]])
Solo_Exposto = data.frame(Classe = "Solo Exposto", atributos[[4]])
Vegetacao = data.frame(Classe = "Vegetação", atributos[[5]])
Vegetacao_Rasteira = data.frame(Classe = "Vegetação Rasteira", atributos[[6]])

amostras.final = rbind(Agricultura, Agua, Area_Urbana, Solo_Exposto, Vegetacao, Vegetacao_Rasteira)

write.csv(amostras.final, 'arquivos_gerados/amostras_extraidas.csv')

#_______________________________________________________________________________

#carregando amostras
amostras = read.table("arquivos_gerados/amostras_extraidas.csv", header = T, sep = ",")[-1]

#calcular espectro de reflectancia medio para cada classe
agrupado = group_by(amostras, Classe)
print(agrupado)

media_ref = summarise_each(agrupado, mean)
print(media_ref)

#calcular a transposta
refs = t(media_ref[,2:10])

cores = c("yellow", "blue", "pink", "brown", "darkgreen", "green")

wavelength = c(490, 560, 660, 705, 740, 770, 850, 1600, 2200)

matplot(x = wavelength, y = refs, type = 'l', lwd = 2, lty = 1,
        xlab = "Comprimento de onda [nm]", ylab = "Reflectância x 10000",
        col = cores, ylim = c(0,8000))
legend('top', legend = media_ref$Classe, lty = 1, col = cores, ncol = 3)

#figura do boxplot
dados.melt = melt(amostras)

ggplot(data = dados.melt, aes(Classe, value, fill = Classe)) +
  geom_boxplot() + 
  facet_wrap(~variable, scale = 'free') + 
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Tahoma"),
        axis.title = element_text(face="bold", size = 10),
        axis.text.x = element_text(colour = "white", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(size = 1, colour = "black")) +
  theme(plot.margin = unit(c(1,1,1,1), "lines"))

#_______________________________________________________________________________

#carregar arquivo com amostras
dados = read.table("arquivos_gerados/amostras_extraidas.csv", header = T, sep = ",")[-1]
head(dados)

#Checar classe (structure) dos dados

str(dados)

dados$Classe = as.factor(dados$Classe)

#separacao dos dados - treinamento e validacao

set.seed(1234)
amostras_treino = sample.split(dados$Classe, SplitRatio = 0.7)

#criar dataframe de treino e teste

train = dados[amostras_treino,]
valid = dados[amostras_treino == F,]

write.csv(train, 'classificacao/dados_treino.csv')
write.csv(valid, 'classificacao/dados_validacao.csv')

#_______________________________________________________________________________

#carregar os arquivos de treino e validacao
treino = read.csv("classificacao/dados_treino.csv", header = T, sep = ",")[-1]
valid = read.csv("classificacao/dados_validacao.csv", header = T, sep = ",")[-1]

#mudar para factor
treino$Classe = as.factor(treino$Classe)
valid$Classe = as.factor(valid$Classe)

#random forest
set.seed(1234)
RF = randomForest(Classe~., data = treino, ntree = 100, mtry = 5, importance = T)
varImpPlot(RF)

importance(RF)

#SVM
set.seed(1234)
SVM = svm(Classe~., kernel = 'polynomial', data = treino) 

pred.RF = predict(RF, valid)
pred.SVM = predict(SVM, valid)


#matriz de confusao
CM.RF = confusionMatrix(data = pred.RF, reference = valid$Classe)
CM.SVM = confusionMatrix(data = pred.SVM, reference = valid$Classe)

print(CM.RF)
print(CM.SVM)

#salvar modelos

saveRDS(object = RF, file = 'classificacao/classificacao_rf.rds')
saveRDS(object = SVM, file = 'classificacao/classificacao_svm.rds')

#_______________________________________________________________________________

#carregar modelos
RF = readRDS("classificacao/classificacao_rf.rds")
SVM = readRDS("classificacao/classificacao_svm.rds")

#carregar TIF com nome das variaveis
img = stack("arquivos_gerados/brasilia_com_indices.tif")

nomes = read.table("arquivos_gerados/nomes_bandas.csv", header = T, sep = ",")

names(img) = nomes[,2]
names(img)

#fazer a predicao para o raster

RF.raster = predict(img, RF)
SVM.raster = predict(img, SVM)

#plot do RF
cores = c("yellow", "blue", "pink", "brown", "darkgreen", "green")
classes = c('Agricultura', 'Agua', 'Área Urbana', 'Solo Exposto', 'Vegetação', 'Vegetação Rasteira')

jpeg(filename = 'classificacao/RF_class.jpeg', width = 15, height = 15, res = 200, units = 'in')
plot(RF.raster, legend = FALSE, col = cores, main = "Classificação Random Forest \n Brasília-DF",
     cex.axis = 1.5, cex.main = 1.5)

legend('topleft', legend = classes, fill = cores, border = FALSE, cex = 2)

addscalebar()
addnortharrow(cols = c("black", "black"), scale = 0.755)

dev.off()

#salvar TIF

writeRaster(RF.raster, 'arquivos_gerados/random_forest_class.tif')
writeRaster(SVM.raster, 'arquivos_gerados/SVM_class.tif')



