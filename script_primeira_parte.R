setwd("C:/Users/Utilizador/Desktop/Universidade/Bioinformática 1º ano/2º Semestre/Extração de Conhecimento de Dados Biológicos/Enunciado do trabalho 1")
getwd()

#instalar o pacote BiocManager caso ele não esteja já instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#instalar o pacote TCGAbiolinks caso ele não esteja já instalado
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

BiocManager::install("SummarizedExperiment", dependencies = TRUE)


#carregar os pacotes na sessão atual do R
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)


#criar uma consulta ao Genomic Data Commons (GDC) e obter dados de perfilamento 
#transcriptômico do projeto TCGA sobre carcinoma endometrial uterino (UCEC)
projeto <- "TCGA-UCEC"
query <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#baixar os dados do Genomic Data Commons (GDC) com base nas especificações 
#de uma consulta criada anteriormente
GDCdownload(query=query)
#preparar os dados baixados do Genomic Data Commons (GDC) para análise
rna_seq_UCEC  <- GDCprepare(query = query, save = TRUE, save.filename = "mRNA_TCGA-UCEC.rda")

#retorna a classe do objeto (tipo de dados ou a estrutura de dados que o objeto 
#representa, o que por sua vez determina quais funções podem ser aplicadas a ele)
class(rna_seq_UCEC)
#retorna as dimensões do objeto(matriz ou um dataframe;objeto mais complexo)
dim(rna_seq_UCEC)

#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_UCEC
geneExp <- SummarizedExperiment::assay(rna_seq_UCEC)
#para ver os outliers - no nosso caso é zero
pre = TCGAanalyze_Preprocessing(rna_seq_UCEC) #faz correlação e tenta identificar outliers
dim(geneExp)
dim(pre)

#NOTA: a normalização pode ser realizada com o package deseq2 e com a função normalization
#NOTA: o edger pode fazer outro processo de normalização através da função rpkm 
#(tem que se fazer a normalização ou antes ou depois da analise diferencial)

#Ver nomes das colunas
colnames(rna_seq_UCEC)
#nomes das linhas
names(rna_seq_UCEC)


#atribui a um novo objeto chamado meta_UCEC, os metadados associados ao conjunto 
#de dados rna_seq_UCEC
meta_UCEC = colData(rna_seq_UCEC)
#retorna as dimensões dos metadados 
dim(meta_UCEC)
#Ver nomes das colunas
colnames(meta_UCEC)
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
#Análise exploratória
#Pré-processamento e filtragem
#retirar as colunas dos metadados onde havia mais de 60 elementos como: 
#“not/Not reported/Reported” e/ou “NA”

cols_with_not_reported <- which(sapply(meta_UCEC,function(x) sum(x == "not reported", na.rm = TRUE)) > 60)
cols_with_Not_Reported <- which(sapply(meta_UCEC,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 60)
cols_with_NA <- which(sapply(meta_UCEC, function(x) sum(is.na(x))) > 60)
# remover as colunas baseadas nos critérios específicos de cima
metadata_matriz_clean <- meta_UCEC[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)] 
dim(metadata_matriz_clean)


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

