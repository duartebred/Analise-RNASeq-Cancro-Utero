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


# Realização de uma consulta ao Genomic Data Commons (GDC) para obter os dados de expressão referentes ao projeto em estudo
projeto <- "TCGA-UCEC"
query_TCGA_UCEC <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")


# Download e preparação dos dados com base na query
GDCdownload(query=query_TCGA_UCEC)
rna_seq_UCEC  <- GDCprepare(query = query_TCGA_UCEC, save = TRUE, save.filename = "mRNA_TCGA-UCEC.rda")


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
geneExp = SummarizedExperiment::assay(rna_seq_UCEC)


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
# vital_status
tabela_vital = prop.table(table(amostras_meta_reduzido$vital_status, useNA = "ifany")) *100
barplot(tabela_vital, names.arg= c("Alive","Dead","NAs"), col= "lightblue", ylab="Percentagem",
        main="Distribuição do estado vital dos Pacientes")


# primary_diagnosis
tabela_primary = prop.table(table(amostras_meta_reduzido$primary_diagnosis, useNA = "ifany")) * 100
outros = sum(tabela_primary[c(1, 2, 3, 5, 6, 8, 9)])
primary_comprimido = matrix(c(outros, tabela_primary[4], tabela_primary[7]), ncol = 3,
    dimnames = list(NULL, c("Outros", "Endometrioid adenocarcinoma, NOS", "Serous cystadenocarcinoma, NOS")))
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


# teste de hipótese
#idade dos pacientes em função do estádio FIGO
boxplot(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage, horizontal = T,
        xlab="Idade do paciente", ylab="Estádio FIGO", main= "Idade do Paciente por Estádio FIGO",
        col=c("lightcoral","indianred","tomato",'red')) #verificação visual dos vários grupos


shapiro.test(metados_sem_nas$age_at_index) # teste à normalidade dos dados
qqnorm(metados_sem_nas$age_at_index) # visualização da normalidade através do gráfico qqplot
qqline(metados_sem_nas$age_at_index)

aov_FIGO = aov(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage); summary(anova_one_way)
boxplot(aov_FIGO$residuals, horizontal = T ) # verificação da homogeneidade das variâncias
t.test(aov_FIGO$residuals, mu=0) # verificação da homogeneidade das variâncias

kruskal.test(metados_sem_nas$age_at_index~metados_sem_nas$figo_stage) # teste não parametrico


#idade dos pacientes em função do estado vital
boxplot(metados_sem_nas$age_at_index~metados_sem_nas$vital_status, horizontal = T,
        xlab="Idade do paciente", ylab="Estado Vital", main= "Idade do Paciente por Estado Vital",
        col=c("lightgreen","gray"))
# a idade continua a não possuir distribuição normal

aov_estado_vital = aov(metados_sem_nas$age_at_index~metados_sem_nas$vital_status); summary(anova_one_way)
boxplot(aov_estado_vital$residuals, horizontal = T ) # verificação da homogeneidade das variâncias
t.test(aov_estado_vital$residuals, mu=0) # verificação da homogeneidade das variâncias

kruskal.test(metados_sem_nas$age_at_index~metados_sem_nas$vital_status) # teste não parametrico

# ALTERNATIVA

# contagens por milhao
gene_exp_filtrado = rna_seq_UCEC[,!is.na(rna_seq_UCEC$figo_stage)]
myCPM = cpm(assay(gene_exp_filtrado), log = T)
head(myCPM)
# filtra dados para ter apenas os genes com mais de 0.5 em pelo menos 2 amostras
thresh = myCPM > 0.5
keep = rowSums(thresh) >= 2
counts.keep = assay(gene_exp_filtrado)[keep,]
summary(keep)
dim(counts.keep)
# cria objeto DGEList
dgeObj = DGEList(counts.keep)

design = model.matrix(~ gene_exp_filtrado$figo_stage)
design

dgeObj = estimateCommonDisp(dgeObj)
dgeObj = estimateGLMTrendedDisp(dgeObj)
dgeObj = estimateTagwiseDisp(dgeObj)

fit <- glmFit(dgeObj, design)
lrt.BvsL <- glmLRT(fit, coef=2)
topTags(lrt.BvsL)


results <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results
dim(results)
summary(de <- decideTestsDGE(lrt.BvsL))


detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)

signif <- -log10(results$FDR)
plot(results$logFC,signif,pch=16)
points(results[detags,"logFC"],-
           log10(results[detags,"FDR"]),pch=16,col="red")


var_genes <- apply(myCPM, 1, var)
select_var <- names(sort(var_genes,
                           decreasing=TRUE))[1:500]
highly_variable_lcpm <- myCPM[select_var,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <-  c("purple","orange")[gene_exp_filtrado$figo_stage]
heatmap.2(highly_variable_lcpm,
            col=rev(morecols(50)),
            trace="column",
            main="Top 500 most variable genes
across samples",
            ColSideColors=col.cell,scale="row")



# filtragem dos dados da contagem de expressão de genes
# filtrar as amostras dos dados de expressão acordo com o estádio FIGO 
gene_exp_filtrado = rna_seq_UCEC[,!is.na(rna_seq_UCEC$vital_status)] 
gene_exp_filtrado$vital_status = factor(gene_exp_filtrado$vital_status)
gene_exp_filtrado = rna_seq_UCEC[,!is.na(rna_seq_UCEC$figo_stage)] 
gene_exp_filtrado$figo_stage = gsub(".*\\b(Stage [VI]+).*", "\\1", gene_exp_filtrado$figo_stage)
gene_exp_filtrado$figo_stage = gsub(" ", "_", gene_exp_filtrado$figo_stage)
gene_exp_filtrado$figo_stage = factor(gene_exp_filtrado$figo_stage) #  transformar o metadado para fator
ddsSE = DESeqDataSet(gene_exp_filtrado, design = ~ figo_stage)
ddsSE = DESeqDataSet(gene_exp_filtrado, design = ~ vital_status)

# filtrar o número de genes
#avaliar se existe com valores omissos
sum(is.na(gene_exp_filtrado)) #não existem valores omissos

#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# segundo a informação contida na documentação do package deseq2 a filtragem de genes com baixa expressão
# mantendo os genes com mais de 10 counts em pelo menos 3 amostras

# Encontrar genes com mais de 10 contagens e que aparecem em mais de 3 amostras
genes_manter = rowSums(counts(ddsSE) > 30)
genes_manter = genes_manter[genes_manter > 4]

genes_manter <- rowSums(counts(ddsSE) >= 30) >= 5
ddsSE <- ddsSE[genes_manter, ]



# Filtrar o DESeqDataSet
ddsSE = ddsSE[genes_manter, ]
dim(ddsSE)


# Teste de expressão diferencial, este package faz a normalização antes de proceder aos testes
# explicar o método de normalização utilizado pelo DESeq2
# 1. Extrair as contagens
counts <- counts(ddsSE)

# 2. Aplicar a transformação logarítmica (por exemplo, log2)
counts_log <- round(log10(counts + 1))  # Adiciona 1 para evitar log(0)

# 3. Substituir as contagens originais pelos valores transformados
assay(ddsSE) <- counts_log


ddsSE_norm = DESeq(ddsSE)
ddsSE_vst <- vst(ddsSE_norm) #transformação de variância estabilizadora das counts pois uma análise inicial revelou elevado número de low reads


resultados = results(ddsSE_norm)
summary(resultados) # sumario dos resultados do teste de expressão diferencial
sum(resultados$padj < 0.1, na.rm=TRUE) # número total de genes diferencialmente expressos


DESeq2::plotMA(resultados, main="DESeq2") # visualização gráfica dos resultados, pontos azuis genes DE


# Análise dos resultados da expressão diferencial


# análise individual do gene mais diferencialmente expresso
plotCounts(ddsSE_norm, gene=which.min(resultados$padj), intgroup="figo_stage", pch = 19, col=1:4)
#NOTA: a normalização pode ser realizada com o package deseq2 e com a função normalization
#NOTA: o edger pode fazer outro processo de normalização através da função rpkm 
#(tem que se fazer a normalização ou antes ou depois da analise diferencial)
