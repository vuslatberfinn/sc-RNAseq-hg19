
# Paketlerin yuklenmesi

 install.packages("Seurat")
 install.packages("tidyverse")
 #SeuratDisk paketi github kodu

library(Seurat)     #Seurat paketi aktif edilir.
library(Seuratdisk)    #SeuratDisk paketi aktif edilir
library (tidyverse)
# RDS format
rds_obj <- readRDS("ependymal_cells.rds")

#10X CellRanger .HDF5 format
#.mtx file
#.loom file
#.h5ad format

#### Veri Setinin Yuklenmesi

pbmc.data <- Read10X(data.dir= "/Users/vuslatberfinsakizci/Desktop/filtered_gene_bc_matrices/hg19")


pbmc<- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells =3 , min.features= 200)
pbmc
# CreateSeuratObject: Seurat paketi ile bir Seurat nesnesi (örnekte pbmc) oluşturur. 
#Bu nesne, scRNA-seq analizleri boyunca verileri taşıyan temel veri yapısıdır.
#mincells:En az 3 hücrede ifade edilen genleri tutar. Bu, yalnızca belirli bir hücre sayısında bulunan genleri analizde kullanarak veri boyutunu küçültmek ve gürültüyü azaltmak için yapılır.
#min.features: Her hücrede en az 200 genin ifade edilmiş olması şartını koyar
#Bu, yeterince bilgi barındırmayan (örneğin, düşük kaliteli veya ölü hücrelerden gelen) hücrelerin analizden çıkarılmasına yardımcı olur.
View(pbmc@meta.data)

##  2-) Standart On-Islem Is akisi

###### QC (Kalite Kontrol) ve Hucre Filtreleme

#Her hücrede tespit edilen benzersiz genlerin sayısı.
# 1)Düşük kaliteli hücreler veya boş damlacıklar genellikle çok az gene sahip olacaktır
# 2)Hücre çiftleri veya çokluları anormal derecede yüksek gen sayısı gösterebilir
# Benzer şekilde, bir hücre içinde tespit edilen moleküllerin toplam sayısı (benzersiz genlerle güçlü bir şekilde ilişkilidir)
# Mitokondriyal genoma eşlenen okumaların yüzdesi
# 1)Düşük kaliteli/ölmekte olan hücreler genellikle kapsamlı mitokondriyal kontaminasyon sergiler
# 2)PercentageFeatureSet()Sayımların bir dizi özellikten kaynaklanan yüzdesini hesaplayan işlevle mitokondriyal QC ölçümlerini hesaplıyoruz
# 3)Tüm genlerin kümesini MT-mitokondriyal genler kümesi olarak kullanarak başlıyoruz

pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

View(pbmc@meta.data)
head(pbmc@meta.data, 5)

#orig.ident: Ornegin kimligi veya ornek adi
#nCount_RNA: Hucre basina toplam RNA okuma sayisi
#nFeature_RNA: Hucre basina tespit edilen gen sayisi
#percent.mt: Hucre basina mitokondriyal RNA'nin toplam RNA'ya orani


#QC metriginin gorsellestirilmesi
VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

#FeatureScatter genellikle ozellik-ozellik iliskilerini gorsellelstirmek icin kullanilir 
#plot1 <- FeatureScatter()

#Kaliteli veriler ile devam etmek icin filtreleme islemi gerceklestirelim

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
#nFeature_RNA < 2500 altinda olan hucreleri secer
#nFeature_RNA > 200 tespit edilen gen sayisi 200 un uzerinde olan hucreleri secer
#percent.mt<5: Mitokondriyal RNA yuzdesi %5 in altinda olan hucreleri secer 
#Bu tur filtreleme genellikle verilerin kalitesini arttirmak icin yapilir

#Filtreleme sonrasi tekrar farki anlamak icin gorsellestirme yapilir
VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"))



### 2) Verilerin Normallestirilmesi
#Istenmeyen hucreleri veri kumesinden kaldirdiktan sonraki adim,verileri normallestirmektir.Varsayilan olarak,her hucre icin
#ozellik ifadesi olcumlerini toplam ifadeye gore normallestiren, bunu bir olcek faktoruyle (varsayilan olarak 10.000)
#carpan ve sonucu logaritmik donusturen kuresel olcekleme normallestirme yontemi "LogNormalize" kullaniriz

pbmc <- NormalizeData(pbmc)

## Yuksek Degisken ozelliklerin Tanimlanmasi (Feature Selection)

#Daha sonra veri setinde yüksek hücre-hücre varyasyonu gösteren özelliklerin bir alt kümesini hesaplıyoruz (yani, bazı hücrelerde
#yüksek oranda ifade edilirken, diğerlerinde düşük oranda ifade ediliyorlar). Biz ve diğerleri, bu genlere akış aşağı analizde
#odaklanmanın tek hücreli veri setlerindeki biyolojik sinyali vurgulamaya yardımcı olduğunu bulduk.

pbmc<- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 2000)
pbmc

#En yuksek 10 Degisken geni secme 
top10 <- head(VariableFeatures(pbmc), 10)
top10

#Etiketli ve Etiketsiz Degisken Ozellikleri Gorsellestirme
plot1<- VariableFeaturePlot(pbmc). #siyah kisimlar herhangi bir ekspresyon seviyesi yaratmadiklari icin sifir kabul olur bu yuzden onemsiz 
plot2<- LabelPoints(plot =plot1,points = top10, repel = TRUE) # repel : ust uste cakismamasi icin kullandik

plot1+plot2



  ####Verilerin ölçeklenmesi (visualization of data)

#Sonra, PCA gibi boyutsal indirgeme tekniklerinden önce standart bir ön işleme adımı olan doğrusal bir dönüşüm ('ölçekleme') uygularız. Fonksiyon ScaleData():
#Her genin ifadesini değiştirir, böylece hücreler arası ortalama ifade 0 olur
#Her genin ifadesini ölçeklendirir, böylece hücreler arası varyans 1 olur
#Bu adım, aşağı akış analizlerinde eşit ağırlık sağlar, böylece yüksek oranda ifade edilen genler baskın olmaz
#Bunun sonuçları şurada saklanır:pbmc[["RNA"]]$scale.data
#Varsayılan olarak yalnızca değişken özellikler ölçeklenir.
#featuresEk özellikleri ölçeklendirmek için argümanı belirtebilirsiniz

#Olcekleme ozelliklerin birbirine benzer bir olcekte olmasini saglayarak,modelleme sureclerinde tutarliligi arttirir.Ozellikler(genler)farklli olcekte
#olabilir.Ornegin bazi genlerin ekspresyon seviyeleri digerlerinden cok daha yuksek olabilir.Bu durum bazi genlerin analizin sonuclarini gereksiz
#yere etkilenmesine neden olabilir.Olcekleme tum genlerin esit agirlikta degerlendirilmesini saglar

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)



## Perform Linear Dimensional Reduction (Dogrusal Boyut Indirgeme)

#Daha sonra ölçeklenmiş veriler üzerinde PCA gerçekleştiririz. Varsayılan olarak, yalnızca önceden belirlenmiş değişken özellikleri
#girdi olarak kullanılır, ancak farklı bir alt küme seçmek isterseniz argüman kullanılarak tanımlanabilir 
#features(özel bir özellik alt kümesi kullanmak isterseniz, bunları first'e ilettiğinizden emin olun ScaleData).

#İlk temel bileşenler için Seurat, veri kümesindeki tek hücreler arasında korelasyon (veya anti-korelasyon) gösteren 
#gen modüllerini temsil eden en pozitif ve negatif yüklemelere sahip genlerin bir listesini çıkarır.


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend(). # ikisinide gorsellestirmek icin kulalandik

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
####
ElbowPlot(pbmc)






### 7) Cluster the Cells

#Hucreleri kumeleme,scRNA-seq veri setindeki hucreleri gen ekspresyonlari veya diger ozelliklerine gore benzerliklerine
#dayanarak gruplandirma surecidir.Bu,veri setindeki farkli hucre turlerini veya durumlarini tanimlamaniza yardimci olur

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5). #resolution degerlerin degisebilir
pbmc <- FindClusters(pbmc, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1))
View(pbmc@meta.data)

#Kumeler Idents() fonksiyonu ile bulunabilir
head(Idents(pbmc),5)

## Run non-linear dimensional reduction (UMAP/tSNE)
#Kumeleme sonuclari genellikle t-SNE veya UMAP ile gorsellestirlir
#Her bir kume farkli renklerle gosterilir,bu da hucrelerin benzerliklerini gore gruplandirilmasini saglar

pbmc <- RunUMAP (pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims= 1:10)


DimPlot(pbmc, reduction = "umap", label = T) #resolution degerlerini degistirerek grafik yaparak sonuclara gore karar veririz
DimPlot(pbmc, reduction = "tsne")
a1<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.0.2")
a2<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.0.4")
a3<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5")
a4<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.0.6")
a5<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.0.8")
a6<-DimPlot(pbmc, reduction = "umap", label = T, group.by = "RNA_snn_res.1")
a1+a2+a3+a4+a5+a6


#Finding differentially expresse features (cluster biomarkers)

#Seurat, diferansiyel ifade (DE) yoluyla kümeleri tanımlayan işaretleyicileri bulmanıza yardımcı olabilir. Varsayılan olarak, tek bir
#kümenin (içinde belirtilen ident.1) pozitif ve negatif işaretleyicilerini, diğer tüm hücrelerle karşılaştırarak belirler.
#FindAllMarkers()bu işlemi tüm kümeler için otomatikleştirir, ancak küme gruplarını birbirlerine veya tüm hücrelere karşı da test edebilirsiniz.

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE) #Bu kodun calismasi uzun suruyor
                                                      #Seurat yeni bir degisken olusturmani istiyor 
                                                    

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
 
#Gorsellestirme

VlnPlot(pbmc, feature= c("MS4A1","CD79A")) #MANUEL sekilde yapiyoruz
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))

## Asigning cell type identity to clusters (kumelere hucre tipi )


new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


FeaturePlot(pbmc,features= c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"), label = T, pt.size = 0.5)
############################################################################################################################################
#################################################################################################################


