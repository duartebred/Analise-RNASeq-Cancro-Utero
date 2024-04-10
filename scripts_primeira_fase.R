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


# Realização de uma consulta ao Genomic Data Commons (GDC) para obter os dados de expressão referentes ao projeto em estudo
projeto <- "TCGA-UCEC"
query_TCGA_UCEC <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")


# Download e preparação dos dados com base na query
GDCdownload(query=query_TCGA_UCEC, method = "api")
rna_seq_UCEC  <- GDCprepare(query = query_TCGA_UCEC, save = FALSE)


#loading dos dados a partir dos ficheiros criados no GDCprepare
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
ddsSE = DESeqDataSet(gene_exp_filtrado, design = ~ vital_status)
dim(ddsSE)


genes_manter = rowSums(counts(ddsSE) >= 15) >= 3
ddsSE = ddsSE[genes_manter, ]
dim(ddsSE)

ddsSE_norm = DESeq(ddsSE)
resultados = results(ddsSE_norm, alpha = 0.05)

summary(resultados) # sumario dos resultados do teste de expressão diferencial
sum(resultados$padj < 0.05, na.rm=TRUE) # número total de genes diferencialmente expressos


DESeq2::plotMA(resultados, main="DESeq2") # visualização gráfica dos resultados, pontos azuis genes DE

plotCounts(ddsSE_norm, gene=10, intgroup="vital_status", pch = 19)

# heatmap
vsd <- varianceStabilizingTransformation(ddsSE_norm, blind = FALSE)
resOrdered[1,]
resOrdered = resultados[order(resultados$padj),]
select = rownames(head(resOrdered,20))
vsd.counts = assay(vsd)[select,]
df = as.data.frame(colData(rna_seq_UCEC)[,"vital_status"])
pheatmap(vsd.counts, cluster_rows=TRUE,show_colnames = F)


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