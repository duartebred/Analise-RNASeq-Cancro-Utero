#instalar o pacote BiocManager caso ele não esteja já instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment", dependencies = TRUE)


#carregar os pacotes na sessão atual do R
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)


# Realização de uma consulta ao Genomic Data Commons (GDC) para obter os dados de expressão referentes ao projeto em estudo
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
geneExp <- SummarizedExperiment::assay(rna_seq_UCEC)


# Metadados das amostras
## Análise exploratória
not_reported = sapply(amostras_metadados,function(x) sum(x == "not reported", na.rm = TRUE))
Not_Reported = sapply(amostras_metadados,function(x) sum(x == "Not Reported", na.rm = TRUE))
nas = sapply(amostras_metadados,function(x) sum(is.na(x)))
par(mfrow = c(1, 3))
hist(not_reported[not_reported!=0],ylab="Número de colunas")
hist(Not_Reported[Not_Reported!=0],ylab="Número de colunas")
hist(nas[nas!=0],ylab="Número de colunas")

# remoção de colunas com mais de 10% de valores omissos
rm_not_reported = which(sapply(amostras_metadados,function(x) sum(x == "not reported", na.rm = TRUE)) > 60)
rm_Not_Reported = which(sapply(amostras_metadados,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 60)
rm_nas = which(sapply(amostras_metadados, function(x) sum(is.na(x))) > 60)
amostras_meta = amostras_metadados[, -c(rm_not_reported, rm_Not_Reported, rm_nas)]


colnames(metadata_matriz_clean)

table(as.data.frame(metadata_matriz_clean$vital_status))
summary(as.data.frame(metadata_matriz_clean$age_at_diagnosis))
pie(table(as.data.frame(metadata_matriz_clean$primary_diagnosis)))
summary(as.data.frame(metadata_matriz_clean$age_at_index))
boxplot(as.data.frame(metadata_matriz_clean$age_at_index),horizontal=T)
hist(metadata_matriz_clean$age_at_index)
table(as.data.frame(metadata_matriz_clean$figo_stage))
summary(as.data.frame(metadata_matriz_clean$figo_stage))
sum(is.na(metadata_matriz_clean$figo_stage))
dim(metadata_matriz_clean)







# criar vetor logica em que diz que linhas é que possuem grade diferente de na
coluna=!is.na(meta_UCEC$paper_tumor_grade) 


# seleciona apenas as linhas que tem grade
meta=meta_UCEC[coluna,c("paper_tumor_grade","paper_age","disease_type")] 
meta$disease_type
table(meta_UCEC[rownames(meta),"paper_tumor_grade"]) # confirmar que os pacientes selecionados tem um valor de grade atribuido na tabela de metadados original
head(meta)

summary(meta$paper_age)
boxplot(meta$paper_age,horizontal=T)

anova=aov(meta$paper_age~meta$paper_tumor_grade)
summary(anova)

TukeyHSD(anova)

boxplot(meta$paper_age~meta$paper_tumor_grade,horizontal=T)
#verificaçáo da dimensão do dataframe meta após selação de apenas as linhas com informação relativa ao grade
sum(table(meta_UCEC$paper_tumor_grade))
dim(meta)


#criação do data frame de expressão apenas para os pacientes com informação relativa ao grade
exp_grade=geneExp[,rownames(meta)]
dim(exp_grade)



row.names(meta_UCEC[,"paper_tumor_grade"])
select=meta_UCEC$paper_tumor_grade
row_names=rownames(select)

#nomes das linhas
names(meta_UCEC)
#extrair componentes de um objeto por nome (através de colunas)
meta_UCEC$patient
meta_UCEC$paper_vital_status


##################################
#Processamento de Dados Clínicos com TCGAbiolinks
## dados clinicos

#buscar dados clínicos do projeto TCGA-UCEC (carcinoma endometrial uterino), 
#especificamente suplementos clínicos, através de GDCquery
query_clin <- GDCquery(project = "TCGA-UCEC", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement", 
                       data.format = "BCR Biotab")
#Baixa e prepara os dados clínicos
GDCdownload(query = query_clin)
clinical.UCEC <- GDCprepare(query = query_clin, save = TRUE, save.filename = "clinical_data_UCEC.rda")
#listar os nomes dos componentes do objeto
names(clinical.UCEC)

#exibir as primeiras linhas da coluna ou lista clinical_drug_ucec contida dentro do objeto clinical.UCEC
head (clinical.UCEC$clinical_drug_ucec)
#converte os dados contidos em clinical.UCEC$clinical_patient_ucec para um 
#dataframe e atribui esse novo dataframe à variável df
df = as.data.frame(clinical.UCEC$clinical_patient_ucec)
View(df)


#######################################





# Extrair apenas a palavra "Stage" e os números romanos
figo_stage <- gsub(".*\\b(Stage [VI]+).*", "\\1", metadata_matriz_clean$figo_stage)


# Visualizar o resultado
print(metadata_matriz_clean$figo_stage)


#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_UCEC
geneExp <- SummarizedExperiment::assay(rna_seq_UCEC)







##########

#Análise de Expressão Diferencial com DESeq2
library(DESeq2)  #Nota: pacote DESeq2, uma ferramenta para análise de expressão diferencial de dados de contagem de sequenciamento de RNA (RNA-Seq)

#Filtra os dados de RNA para incluir apenas amostras não nulo
data_de <- rna_seq_UCEC[,!is.na(rna_seq_UCEC$paper_vital_status)]

#Cria um objeto DESeqDataSet para análise, especificando um design experimental 
#que compara o status de IDH
ddsSE <- DESeqDataSet(data_de, design = ~ paper_vital_status)

#Filtragem de Genes: Remove genes com contagens baixas (menos de 10) 
#para melhorar a confiabilidade da análise de expressão diferencial.
keep <- rowSums(counts(ddsSE)) >= 10

#Executa a análise de expressão diferencial com a função DESeq
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)
#comparação "WT vs Mutant" para o status do IDH, e converte os resultados para um dataframe
res <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant")
dea <- as.data.frame(res)
#resume os resultados para obter uma visão geral dos achados estatísticos, 
#como o número de genes significativamente diferencialmente expressos
summary(res)



#para ver os outliers - no nosso caso é zero
pre = TCGAanalyze_Preprocessing(rna_seq_UCEC) #faz correlação e tenta identificar outliers
dim(geneExp)
dim(pre)

#NOTA: a normalização pode ser realizada com o package deseq2 e com a função normalization
#NOTA: o edger pode fazer outro processo de normalização através da função rpkm 
#(tem que se fazer a normalização ou antes ou depois da analise diferencial)
