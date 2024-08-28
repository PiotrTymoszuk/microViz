# Pre-processing of the Metastatic Breast Cancer Project gene expression data
# which will be used in examples for the microViz package

# tools -------

  library(tidyverse)
  library(rlang)
  library(stringi)

# reading the expression data ------

  expression <-
    read_tsv('./inst/examples/brca_mbcproject_2022/data_mrna_seq_v2_rsem.txt')

  ## selection of genes with the symbol assigned

  expression <- expression %>%
    filter(!is.na(Hugo_Symbol),
           !duplicated(Hugo_Symbol)) %>%
    select(-Entrez_Gene_Id) %>%
    column_to_rownames('Hugo_Symbol') %>%
    t %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

# clinical data ---------

  clinics <-
    read_tsv('./inst/examples/brca_mbcproject_2022/data_clinical_sample.txt',
             skip = 4)

  ## extracting the sample ID, time point, histology and ER status

  clinics <- clinics %>%
    transmute(sample_id = SAMPLE_ID,
              patient_id = PATIENT_ID,
              timepoint = stri_extract(SAMPLE_TIMEPOINT, regex = 'T\\d{1}$'),
              timepoint = factor(timepoint,
                                 c('T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T7')),
              metastasis = factor(CALC_MET_SETTING,
                                  c('NO_METASTATIC_DISEASE_PRESENT',
                                    'METASTATIC_DISEASE_PRESENT')),
              histology = factor(BX_HISTOLOGY,
                                 c('ADENOCARCINOMA',
                                   'IDC',
                                   'MIXED_IDLC',
                                   'ILC',
                                   'CARCINOMA',
                                   'DCIS')),
              er_status = factor(BX_ER, c('POSITIVE', 'NEGATIVE')))

# an example dataset -------

  brca <- right_join(clinics, expression, by = 'sample_id')

  save(brca, file = './data/brca.RData')

# END ------
