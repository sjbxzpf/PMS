load("AIBL_fea.RData")
load("brain_fea.RData")
load("go_BP.RData")
library(BioM2)
library(mlr3verse)#查看可用学习器mlr_learners
library(caret)
library(parallel)
set.seed(666)
load("AIBL_beta.RData")
AIBL_PATH <- BioM2(TrainData = AIBL_beta,
                  TestData = NULL,
                  pathlistDB = go_BP,
                  FeatureAnno = AIBL_fea,
                  nfolds = 10,
                  classifier = 'liblinear',
                  predMode = "probability",
                  PathwaySizeUp = 200,
                  PathwaySizeDown = 20,
                  MinfeatureNum_pathways = 10,
                  Add_UnMapped = T,
                  Unmapped_num = 30,
                  Add_FeartureSelection_Method = "cor",
                  Inner_CV = T,
                  inner_folds = 5,
                  Stage1_FeartureSelection_Method = 'cor',
                  cutoff = 0,
                  Stage2_FeartureSelection_Method = "cor",
                  cutoff2 = 0,
                  classifier2 = NULL,
                  target = 'pathways',
                  p.adjust.method = "fdr",
                  cores = 5
   )
save(AIBL_PATH, file = "AIBL_PATH.RData")
load("froCortex_beta.RData")
froCortex_PATH <- BioM2(TrainData = froCortex_beta,
                   TestData = NULL,
                   pathlistDB = go_BP,
                   FeatureAnno = brain_fea,
                   nfolds = 5,
                   classifier = 'liblinear',
                   predMode = "probability",
                   PathwaySizeUp = 400,
                   PathwaySizeDown = 10,
                   MinfeatureNum_pathways = 10,
                   Add_UnMapped = T,
                   Unmapped_num = 30,
                   Add_FeartureSelection_Method = "cor",
                   Inner_CV = T,
                   inner_folds = 3,
                   Stage1_FeartureSelection_Method = 'cor',
                   cutoff = 0,
                   Stage2_FeartureSelection_Method = 'RemoveHighcor',
                   cutoff2 = 0.8,
                   classifier2 = 'liblinear',
                   target = 'pathways',
                   p.adjust.method = "fdr",
                   cores = 10
    )
save(froCortex_PATH, file = "froCortex_PATH2.RData")
load("Neuron_beta.RData")
Neuron_PATH <- BioM2(TrainData = Neuron_beta,
                  TestData = NULL,
                  pathlistDB = go_BP,
                  FeatureAnno = brain_fea,
                  nfolds = 5,
                  classifier = 'liblinear',
                  predMode = "probability",
                  PathwaySizeUp = 200,
                  PathwaySizeDown = 20,
                  MinfeatureNum_pathways = 10,
                  Add_UnMapped = T,
                  Unmapped_num = 30,
                  Add_FeartureSelection_Method = "cor",
                  Inner_CV = T,
                  inner_folds = 2,
                  Stage1_FeartureSelection_Method = 'cor',
                  cutoff = 0,
                  Stage2_FeartureSelection_Method = 'cor',
                  cutoff2 = 0,
                  classifier2 = 'liblinear',
                  target = 'pathways',
                  p.adjust.method = "fdr",
                  cores = 10
   )
save(Neuron_PATH, file = "Neuron_PATH.RData")
load("Glia_beta.RData")
Glia_PATH <- BioM2(TrainData = Glia_beta,
                  TestData = NULL,
                  pathlistDB = go_BP,
                  FeatureAnno = brain_fea,
                  nfolds = 5,
                  classifier = 'liblinear',
                  predMode = "probability",
                  PathwaySizeUp = 200,
                  PathwaySizeDown = 20,
                  MinfeatureNum_pathways = 10,
                  Add_UnMapped = T,
                  Unmapped_num = 30,
                  Add_FeartureSelection_Method = "cor",
                  Inner_CV = T,
                  inner_folds = 2,
                  Stage1_FeartureSelection_Method = 'cor',
                  cutoff = 0,
                  Stage2_FeartureSelection_Method = 'cor',
                  cutoff2 = 0,
                  classifier2 = 'liblinear',
                  target = 'pathways',
                  p.adjust.method = "fdr",
                  cores = 10
   )
save(Glia_PATH, file = "Glia_PATH.RData")
ADNI_ROC <- BioM2(TrainData = ADNI_beta,
                  TestData = NULL,
                  pathlistDB = go_BP,
                  FeatureAnno = ADNI_fea,
                  nfolds = 10,
                  classifier = 'liblinear',
                  predMode = "probability",
                  PathwaySizeUp = 200,
                  PathwaySizeDown = 20,
                  MinfeatureNum_pathways = 10,
                  Add_UnMapped = T,
                  Unmapped_num = 100,
                  Add_FeartureSelection_Method = "cor",
                  Inner_CV = T,
                  inner_folds = 5,
                  Stage1_FeartureSelection_Method = 'cor',
                  cutoff = 0,
                  Stage2_FeartureSelection_Method = "RemoveHighcor",
                  cutoff2 = 0.9,
                  classifier2 = NULL,
                  target = 'predict',
                  cores = 10
)
save(ADNI_ROC, file = "ADNI_ROC.RData")