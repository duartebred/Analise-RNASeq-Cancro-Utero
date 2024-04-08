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
geneExp <- SummarizedExperiment::assay(rna_seq_UCEC)


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
        main="Distribuição de Vitalidade dos Pacientes")


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
    'yellow','lightpink','orchid'), main= "Distribuição de Estágios de cancro Ginecológico (FIGO)")


# age_at_index
idade_pacientes = amostras_meta_reduzido$age_at_index
summary(idade_pacientes)
par(mfrow = c(2, 1))
boxplot(idade_pacientes,horizontal=T, col = "purple", main='Histograma da idade dos pacientes')
hist(idade_pacientes, xlab = 'Idade dos Pacientes', ylab = 'Frequência', main = 'Histograma da Idade dos Pacientes',
     col='lightblue')


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
