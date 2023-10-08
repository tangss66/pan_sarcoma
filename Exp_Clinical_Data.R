clinical_data <- clinfo1[,c("firmiana_T","ID","firmiana_N","firmiana_phos_T","firmiana_phos_N","patient",'gender','age','sarcoma.subtype2',
                            'os.time','os.livingstatus','dfs.time','dfs.livingstatus','ki67.pos.ratio','ccp.group_3')]

# clinical_data <- clinical_data[-which(is.na(clinical_data$histological_type)),]
clinical_data <- merge(clinical_data, NumInfo[,c('ID','histological_type')], by = 'ID')

clinical_data$ccp.group_3 <- paste('C',clinical_data$ccp.group_3, sep = '')
ExpData_pro_t <- ExpData[,clinical_data$firmiana_T]
ExpData_pro_t_id <- ExpDataID[,clinical_data$ID] 

ExpData_pho_t <- ExpPhoSite_T
ExpData_pho_n <- ExpPhoSite_N

ExpData_pro_n <- ExpDataN

ExpData_pro_all <- ExpDataTN


save(clinical_data,ExpData_pro_t,ExpData_pro_t_id,ExpData_pro_n,ExpData_pro_all,
     ExpData_pho_t,ExpData_pho_n,file = 'Exp_Clinical_Data.RData')

# 
# colnames(clinfo1)
# [1] "firmiana_T"               "ID"                       "firmiana_N"               "firmiana_phos_T"         
# [5] "firmiana_phos_N"          "patient"                  "Exp_best_N"               "gender"                  
# [9] "age"                      "diagnosis.conclusion"     "sarcoma.subtype1"         "follow.up"               
# [13] "hospital.num"             "diagnosis.date"           "surgery.date"             "death.cause"             
# [17] "living.status"            "death.date"               "radiotherapy"             "radiotherapy.method"     
# [21] "target.therapy"           "target.drug"              "chemotherapy"             "chemotherapy.drug"       
# [25] "recurrence"               "recurrence.times"         "first.recurrence.date"    "second.recurrence.date"  
# [29] "third.recurrence.date"    "forth.recurrence.date"    "metastasis"               "metastasis.times"        
# [33] "first.metastasis.date"    "first.metastasis.region"  "second.metastasis.date"   "second.metastasis.region"
# [37] "third.metastasis.date"    "third.metastasis.region"  "focus"                    "recurrence.date.hospital"
# [41] "follow.up.time"           "lesion"                   "lesion2"                  "sarcoma.subtype2"        
# [45] "malignancy"               "sarcoma.subtype3"         "sarcoma.subtype4"         "death.cause1"            
# [49] "death.date2"              "os.time"                  "os.livingstatus"          "dfs.endtime"             
# [53] "dfs.time"                 "dfs.livingstatus"         "ki67"                     "ki67.pos.ratio"          
# [57] "ccp.group_2"              "ccp.group_3"              "ccp.group_4"              "ccp.group_5"             
# [61] "ccp.group_6"              "ccp.group_7"              "ccp.group_8"              "ccp.group_9"             
# [65] "ccp.group_10"             "os.livingstatus.3y"       "os.time.3y"  
# 
# 
# colnames(NumInfo)
# [1] "ID"                   "HisNum"               "firmiana_T"           "GeneNum_T"           
# [5] "firmiana_N"           "GeneNum_N"            "firmiana_phos_T"      "phospho_peptideNum_T"
# [9] "firmiana_phos_N"      "phospho_peptideNum_N" "gender"               "age"                 
# [13] "os.time"              "ki67.pos.ratio"       "histological_type"    "lesion"              
# [17] "living.status"        "treatment"            "recurrence"           "metastasis"          
# [21] "subtype"             