#instalar o pacote BiocManager caso ele não esteja já instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment", dependencies = TRUE)
BiocManager::install("DESeq2")
BiocManager::install(c("edgeR"))
BiocManager::install(c("Glimma"))
BiocManager::install(c("gplots"))
BiocManager::install(c("org.Mm.eg.db"))
BiocManager::install(c("gplots"))
BiocManager::install("genefilter")
install.packages("Rtsne")
BiocManager::install(c("fgsea"))


#carregar os pacotes na sessão atual do R
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(fgsea)
library(pheatmap)
library(genefilter)
library(Rtsne)
library(fgsea)
library(GSEABase)


# Realização de uma consulta ao Genomic Data Commons (GDC) para obter os dados de expressão referentes ao projeto em estudo
projeto <- "TCGA-UCEC"
query_TCGA_UCEC <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")


# Download e preparação dos dados com base na query
GDCdownload(query=query_TCGA_UCEC, method = "api")
rna_seq_UCEC  <- GDCprepare(query = query_TCGA_UCEC, save = TRUE, save.filename = "mRNA_TCGA-UCEC.rda")


#loading dos dados a partir dos ficheiros criados no GDCprepare
rna_seq_UCEC= get(load("C:/Users/Asus/Downloads/mRNA_TCGA-UCEC.rda"))
#loading ricardo
rna_seq_UCEC= get(load("C:/Users/ricar/Downloads/mRNA_TCGA-UCEC.rda"))

# analise da estrutura dos dados descarregados
class(rna_seq_UCEC)
dim(rna_seq_UCEC)
colnames(rna_seq_UCEC)
rownames(rna_seq_UCEC)


# análise dos metadados
linhas_metadados=SummarizedExperiment::rowData(rna_seq_UCEC)
amostras_metadados = SummarizedExperiment::colData(rna_seq_UCEC)

class(linhas_metadados)
class(amostras_metadados)

dim(linhas_metadados)
dim(amostras_metadados)

colnames(linhas_metadados)
colnames(amostras_metadados)


# Extração da informação relacionada à contagem da expressão dos genes do objeto rna_seq_UCEC
geneExp = SummarizedExperiment::assay(rna_seq_UCEC, "unstranded")



# Metadados das amostras
# remoção de colunas com mais de 10% de valores omissos
rm_not_reported = which(sapply(amostras_metadados,function(x) sum(x == "not reported", na.rm = TRUE)) > 60)
rm_Not_Reported = which(sapply(amostras_metadados,function(x) sum(x == "Not Reported",na.rm = TRUE)) > 60)
rm_nas = which(sapply(amostras_metadados, function(x) sum(is.na(x))) > 60)
amostras_meta_reduzido = amostras_metadados[, -c(rm_not_reported, rm_Not_Reported, rm_nas)]
dim(amostras_meta_reduzido)


# seleção apenas das colunas de interesse e transformação da coluna figo_stage
amostras_meta_reduzido = amostras_meta_reduzido[,c("vital_status","primary_diagnosis","age_at_index","figo_stage")]
amostras_meta_reduzido$figo_stage = gsub(".*\\b(Stage [VI]+).*", "\\1", amostras_meta_reduzido$figo_stage)


# Análise exploratória dos metadados
#descrição das colunas dos metadados selecionadas
table(amostras_meta_reduzido$vital_status)
table(is.na(amostras_meta_reduzido$vital_status))
table(amostras_meta_reduzido$primary_diagnosis)
table(is.na(amostras_meta_reduzido$primary_diagnosis))
table(amostras_meta_reduzido$age_at_index)
table(is.na(amostras_meta_reduzido$age_at_index))
table(amostras_meta_reduzido$figo_stage)
table(is.na(amostras_meta_reduzido$figo_stage))
      
      
# vital_status
tabela_vital = round(prop.table(table(amostras_meta_reduzido$vital_status, useNA = "ifany")) *100)
barplot(tabela_vital, names.arg= c("Alive","Dead","NAs"), col= "lightblue", ylab="Percentagem",
main="Distribuição do estado vital dos Pacientes")
      
      
# primary_diagnosis
tabela_primary = prop.table(table(amostras_meta_reduzido$primary_diagnosis, useNA = "ifany")) * 100
outros = sum(tabela_primary[c(1, 2, 3, 5, 6, 8, 9)])
primary_comprimido = matrix(c(outros, tabela_primary[4], tabela_primary[7]), ncol = 3, dimnames = list(NULL, c("Outros", "Endometrioid adenocarcinoma, NOS", "Serous cystadenocarcinoma, NOS")))
      labels_primary = c("Outros", "Endometrioid adenocarcinoma, NOS", "Serous cystadenocarcinoma, NOS")
      pie(primary_comprimido, labels = paste(labels_primary, sprintf("%.1f%%", primary_comprimido)),
          col = c('lightgreen', 'lightblue', 'lightpink'), main = "Distribuição de Diagnósticos Primários")
      
      
# figo_stage
tabela_figo = prop.table(table(amostras_meta_reduzido$figo_stage,  useNA = "ifany")) *100
labels_figo = c("Stage I", "Stage II", "Stage III", "Stage IV", "NAs")
pie(tabela_figo, labels = paste(labels_figo, sprintf("%.1f%%", tabela_figo)), col=c('lightblue','lightgreen',
                                                                                          'yellow','lightpink','orchid'), main= "Distribuição de estádios de cancro Ginecológico (FIGO)")
      
      
# age_at_index
idade_pacientes = amostras_meta_reduzido$age_at_index
summary(idade_pacientes)
boxplot(idade_pacientes,horizontal=T, col = "purple", main='Histograma da idade dos pacientes',)
hist(idade_pacientes, xlab = 'Idade dos Pacientes', ylab = 'Frequência', main = 'Distribuição da Idade dos Pacientes',
           col='lightblue')
      
      
# eliminação das linhas que possuem NAs
metados_sem_nas = na.omit(as.data.frame(amostras_meta_reduzido))
dim(metados_sem_nas) # passamos de 589 para 574
      
#**Nota:**para todos os testes estatísticos abaixo realizados, considere-se um valor de prova igual a 0.05.
      
#idade dos pacientes em função do estádio FIGO
boxplot(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage, horizontal = T,
        xlab="Idade do paciente", ylab="Estádio FIGO", main= "Idade do Paciente por Estádio FIGO",
        col=c("lightcoral","indianred","tomato",'red')) #verificação visual dos vários grupos
      
      
shapiro.test(metados_sem_nas$age_at_index) # teste à normalidade dos dados
qqnorm(metados_sem_nas$age_at_index) # visualização da normalidade através do gráfico qqplot
qqline(metados_sem_nas$age_at_index)
      
# Teste de Bartlett para verificar homogeneidade das variâncias
#H0: as variâncias são homogénas
#H1: as variâncias não são homogénas
bartlett.test(metados_sem_nas$age_at_index ~ metados_sem_nas$figo_stage)
      
#H0: Não há diferença estatisticamente significativa na média das idades dos pacientes entre os diferentes estádios FIGO
#H1: Há diferença estatisticamente significativa na média das idades dos pacientes entre os diferentes estádios FIGO
aov_FIGO = aov(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage)
summary(aov_FIGO)
boxplot(aov_FIGO$residuals, horizontal = T ) # verificação da homogeneidade das variâncias
      
t.test(aov_FIGO$residuals, mu=0) # verificação da homogeneidade das variâncias
kruskal.test(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage) # teste não parametrico
      
      
#idade dos pacientes em função do estado vital
# Teste de Bartlett para verificar homogeneidade das variâncias
#H0: as variâncias são homogénas
#H1: as variâncias não são homogénas
bartlett.test(metados_sem_nas$age_at_index ~ metados_sem_nas$vital_status)
boxplot(metados_sem_nas$age_at_index~metados_sem_nas$vital_status, horizontal = T,
        xlab="Idade do paciente", ylab="Estado Vital", main= "Idade do Paciente por Estado Vital",
        col=c("lightgreen","gray"))
      
      
aov_estado_vital = aov(metados_sem_nas$age_at_index~metados_sem_nas$vital_status); summary(aov_estado_vital)
boxplot(aov_estado_vital$residuals, horizontal = T ) # verificação da homogeneidade das variâncias
      
t.test(aov_estado_vital$residuals, mu=0) # verificação da homogeneidade das variâncias
kruskal.test(metados_sem_nas$age_at_index~metados_sem_nas$vital_status) # teste não parametrico


# filtrar as amostras que são de "Endometrioid adenocarcinoma, NOS"
amostras_filtradas = amostras_metadados[!is.na(amostras_metadados$primary_diagnosis),]
amostras_filtradas = amostras_filtradas[amostras_filtradas$primary_diagnosis =="Endometrioid adenocarcinoma, NOS",]
dados_EA = geneExp[,rownames(amostras_filtradas)]
dim(dados_EA)

sum(is.na(amostras_filtradas$vital_status))
amostras_filtradas$vital_status = factor(amostras_filtradas$vital_status)

ddsSE = DESeqDataSetFromMatrix(countData = dados_EA, colData = amostras_filtradas, design = ~vital_status)
dim(ddsSE)


# filtragem de genes com menos de 15 ocorrências em pelo menos 3 amostras
genes_manter = rowSums(counts(ddsSE) >= 20) >= 4
ddsSE = ddsSE[genes_manter, ]
dim(ddsSE)

ddsSE_norm = DESeq(ddsSE)
resultados = results(ddsSE_norm, alpha = 0.05)

summary(resultados) # sumario dos resultados do teste de expressão diferencial
sum(resultados$padj < 0.05, na.rm=TRUE) # número total de genes diferencialmente expressos


# visualização gráfica dos resultados
DESeq2::plotMA(resultados, main="Análise de expressão diferencial DESeq2")
plotCounts(ddsSE_norm, gene=which.min(resultados$padj), intgroup="vital_status", pch = 19)


# heatmap
vsd <- varianceStabilizingTransformation(ddsSE_norm, blind = FALSE)
resOrdered = resultados[order(resultados$padj),]
select = rownames(head(resOrdered,20))
vsd.counts = assay(vsd)[select,]
df = as.data.frame(colData(rna_seq_UCEC)[,"vital_status"])
pheatmap(vsd.counts, cluster_rows=TRUE,show_colnames = F)


# análise expressão diferencial através do edgeR (https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
#Pre processamento, normalização contagens por milhão e filtragem de genes com baixa expressão ou expressão nula

# normalização dos dados das contagens
dados_EA_CPM = cpm(dados_EA)
dim(dados_EA_CPM)

# remoção dos genes que não possuem expressão de 0.5 ou superior em pelo menos 2 amostras 
min_exp = dados_EA_CPM > 0.5
keep_rows =  rowSums(min_exp) >= 2
dados_EA_CPM_filt = dados_EA_CPM[keep_rows,]
dim(dados_EA_CPM_filt)


# transformação das contagens normalizadas
dgeObj = DGEList(dados_EA_CPM_filt)
geneExp_filt_log = cpm(dgeObj, log=TRUE)
boxplot(geneExp_filt_log[,1:100], xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(geneExp_filt_log),col="blue")
title("Boxplots das logCPMs (não normalizado)")


# normalização do tamanho da libraria de cada amostra
dgeObj = calcNormFactors(dgeObj)


# análise de expressão diferencial
vital_status = as.factor(amostras_filtradas$vital_status)
design = model.matrix(~vital_status)


# Para estimar a dispersão comum, as dispersões com tendência e as dispersões por tag em uma única execução
dgeObj = estimateDisp(dgeObj, design = design)


fit = glmFit(dgeObj, design)

lrt = glmLRT(fit, coef=2)
topTags(lrt)


# informação acerca dos genes (subexpressos, não significativo, sobreexpresso)
summary(decideTests(lrt))


# visualização gráfica dos genes diferencialmente expressos
plotMD(lrt)
abline(h=c(-1, 1), col="blue")


# glmTreat é usado para filtrar os genes DE e focar apenas nos genes que são mais biologicamente significativos
# filtrar os genes diferencialmente expressos com maior importância biológica
filtered_results = glmTreat(fit, lfc=log2(1.5))
topTags(filtered_results)


summary(decideTests(filtered_results))
plotMD(filtered_results)
abline(h=c(-1, 1), col="blue")


# análise da expressão diferencial dos genes mais comumente mutados https://tcr.amegroups.org/article/view/46888/html
#PTEN(>77%), PIK3CA (53%), PIK3R1 (37%), CTNNB1 (36%), ARID1A (35%), K-RAS (24%), CTCF (20%)
#RPL22 (12%), TP53 (11%), FGFR2 (11%), and ARID5B (11%).

genes_mutados = c("ENSG00000171862.11", "ENSG00000121879.6", "ENSG00000145675.15", "ENSG00000168036.18",
                  "ENSG00000117713.20", "ENSG00000133703.13", "ENSG00000102974.16", "ENSG00000116251.11",
                  "ENSG00000141510.18", "ENSG00000066468.23", "ENSG00000150347.16")


n_gene_edgeR=c()
for (gene in genes_mutados){
  n_gene_edgeR=c(n_gene_edgeR, which(rownames(filtered_results)==gene))
}
filtered_results[n_gene_edgeR,]$table


n_gene_DESeq2=c()
for (gene in genes_mutados){
  n_gene_DESeq2=c(n_gene_DESeq2, which(rownames(resultados)==gene))
}
resultados[n_gene_DESeq2,]


# analise dos genes do estudo sinalizados como diferencialmente expressos
plotCounts(ddsSE_norm, gene=which(rownames(resultados)=="ENSG00000121879.6"), intgroup="vital_status", pch = 19)
plotCounts(ddsSE_norm, gene=which(rownames(resultados)=="ENSG00000133703.13"), intgroup="vital_status", pch = 19)


# enrequecimento Ricardo

# extração do nome de cada gene
nome_gene=c()
for (gene in rownames(resultados)) {
  nome_gene =c(nome_gene,linhas_metadados[gene,"gene_name"])
}


# adição do nome de cada gene ao dataframe dos resultados da expressão diferencial
resultados["gene_name"] =nome_gene

# salvar o data frame resultados
write.csv(resultados, file = "dge_deseq2.csv", row.names = T)

# carregamento do grupo de genes para a análise de enriquecimento
path = gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Ordena os resultados pela alteração na expressão em ordem decrescente
results_ord = resultados[order(-resultados[,"log2FoldChange"]), ]
# Prepara os rankings para a FGSEA
ranks = results_ord$log2FoldChange
# associa às linhas o nome do gene
names(ranks) <- results_ord$gene_name


fgseaRes = fgsea(path, stats = ranks, minSize = 15, maxSize = 500)
head(fgseaRes[order(padj), ])


#Enriquecimento
get_entrez <- function(x) {
  unlist(strsplit(x, split="[.]+"))[2]
}
#Busca anotações de genes e exibe as primeiras linhas do objeto
ann <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = sapply(rownames(resultados), get_entrez), 
  columns = c("ENTREZID", "SYMBOL", "GENENAME"))
head(ann)


# Combina os resultados da expressão diferencial com as anotações de genes
all_results.annotated <- cbind(as.data.frame(resultados), ann)
head(all_results.annotated)


# Ordena os resultados pela alteração na expressão em ordem decrescente
results.ord <- all_results.annotated[order(-all_results.annotated[,"log2FoldChange"]), ]
# Prepara os rankings para a FGSEA
ranks <- results.ord$log2FoldChange
names(ranks) <- results.ord$ENTREZID


pathways <- gmtPathways("C:/Users/Utilizador/Desktop/Universidade/Bioinformática 1º ano/2º Semestre/Extração de Conhecimento de Dados Biológicos/Enunciado do trabalho 1/h.all.v2023.2.Hs.entrez.gmt")


# Executa a FGSEA
fgseaRes <- fgsea(pathways, 
                  stats = vetor,
                  scoreType = 'std',
                  minSize = 10, 
                  maxSize = 1000)


class(fgseaRes)
dim(fgseaRes) 


# Mostra as primeiras entradas ordenadas pelo p-valor ajustado
head(fgseaRes[order(fgseaRes$padj), ])


# Cria um gráfico de barras dos resultados da FGSEA
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA")

#PCA
# Os PCs com menor variância são descartados para reduzir efetivamente a dimensionalidade dos dados sem perder informações.
#Ela faz isso encontrando novas variáveis, chamadas componentes principais, que explicam a maior parte da variação nos dados 
#originais. Essas novas variáveis são criadas por combinações lineares das variáveis originais.

pcares1 = prcomp(dados_EA_CPM, scale = F)    #já estao normalizados
summary(pcares1)$importance[3, ]

"encontrar o número de componentes principais necessários para explicar pelo menos 95% da variância dos dados."

min(which(summary(pcares1)$importance[3,]>0.95))
pcares1$rotation[, 1:20]

plot(pcares1)
biplot(pcares1)


amostras_filtradas$figo_stage = gsub(".*\\b(Stage [VI]+).*", "\\1", amostras_filtradas$figo_stage) #recuperar transformação inicial
amostras_filtradas$figo_stage <- factor(amostras_filtradas$figo_stage)
cores_estagio <- rainbow(length(levels(amostras_filtradas$figo_stage)))
plot(pcares$x, col = cores_estagio[amostras_filtradas$figo_stage], pch = 19)


amostras_filtradas$vital_status <- factor(amostras_filtradas$vital_status)
plot(pcares$x, col=as.integer(amostras_filtradas$vital_status), pch = 19)

#tSNE

Rtsne(dados_EA_CPM)
dados_EA_CPM_nd = dados_EA_CPM[!duplicated(dados_EA_CPM),]
dim(dados_EA_CPM_nd)
res_tnse = Rtsne(dados_EA_CPM_nd)
plot(res_tnse$Y, col = cores_estagio[amostras_filtradas$figo_stage], pch = 19)

amostras_filtradas$vital_status <- factor(amostras_filtradas$vital_status)
plot(res_tnse$Y, col=as.integer(amostras_filtradas$vital_status), pch = 19)

##Clustering

#Hierárquico 

ddsSE_norm <- DESeq(ddsSE) 
data_rna_UCEC_matrix <- as.matrix(assay(ddsSE_norm))
data_rna_UCEC_transposed <- t(data_rna_UCEC_matrix)

#calculo da matrix euclidiana
tt_mdr = rowttests(t(data_rna_UCEC_matrix)) #teste estatístico para cada gene na matriz de dados de RNA. Determina se há diferenças significativas na expressão gênica entre diferentes condições ou grupos experimentais
rank_mdr = order(tt_mdr$p.value) #ranking crescente dos genes com base nos seus pvalues
#NOTA: Quem chegar aqui, não esquecer de adicionar no rmarkdow a libraria do package genefilter
genes_mdr = rank_mdr[1:30] #seleciona os índices dos 30 menores valores-p
data_rna_UCEC_rank = data_rna_UCEC_matrix[genes_mdr,] #cria uma nova matriz que contém apenas as linhas correspondentes aos 30 genes mais diferencialmente expressos

eucD = dist(data_rna_UCEC_rank)


#por genes 
#(usamos a matriz original - data_rna_LGG_matrix(linhas: genes; colunas: amostras))
tt_mdr_g = rowttests(data_rna_UCEC_matrix)
rank_de_mdr_g = order(tt_mdr_g$p.value)
genes_de_mdr_g = rank_de_mdr_g[1:30]
data_rna_UCEC_rank = data_rna_UCEC_matrix[genes_mdr,]
eucD = dist(data_rna_UCEC_rank)

#complete
cl.hier <- hclust(eucD)
plot(cl.hier,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:complete, distância Euclidiana")

#single
cl.hier2 <- hclust(eucD, method="single")
plot(cl.hier2,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:single, distância Euclidiana")

#average
cl.hier3 <- hclust(eucD, method="average")
plot(cl.hier3,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:average, distância Euclidiana")

#heatmap para os 30 genes
heatmap(data_rna_UCEC_rank, labCol = F)



#por paciente 
#(usamos a matriz transposta - data_rna_LGG_transposed(vice-versa), útil quando o foco é amostras)

tt_mdr = rowttests(t(data_rna_UCEC_matrix))
rank_de_mdr = order(tt_mdr$p.value)
genes_de_mdr = rank_de_mdr[1:30]
data_rna_UCEC_rank = data_rna_UCEC_transposed[genes_mdr,]
eucD = dist(data_rna_UCEC_rank)

#complete
cl.hier4 <- hclust(eucD)
plot(cl.hier4,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:complete, distância Euclidiana")

#single
cl.hier5 <- hclust(eucD, method="single")
plot(cl.hier5,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:single, distância Euclidiana")

#average
cl.hier6 <- hclust(eucD, method="average")
plot(cl.hier6,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:average, distância Euclidiana")

#heatmap para os 30 genes
heatmap(data_rna_UCEC_rank, labCol = F)




#k-means

ofs <- c()
for (k in 2:10) {
  kmeans <- kmeans(t(data_rna_UCEC_matrix), centers = k, nstart = 10)
  ofs <- c(ofs, kmeans$tot.withinss)
}
plot_data <- data.frame(num_clusters = 2:10, wss = ofs)

ggplot(plot_data, aes(x = num_clusters, y = wss)) +
  geom_line() +
  geom_point() +
  labs(x = "Num Clusters", y = "WSS") +
  theme_minimal()

#compara os clusters resultantes com uma variável categórica dos dados  (Qual o K??)
resKmeans <- kmeans(t(data_rna_UCEC_matrix),centers=6)
centroides=resKmeans$cluster
table_result=table(centroides, ddsSE$vital_status)
table_result




