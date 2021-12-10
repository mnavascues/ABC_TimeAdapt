library("reticulate")
library("abcrf")

pandas=import("pandas")

project = "test2"
batch = 1
ref_table_sumstats_file = paste0("results/",project,"/",batch,"/ref_table_sumstats.pkl")
ref_table_latent_variables_file = paste0("results/",project,"/",batch,"/ref_table_latent_variables.pkl")


ref_table_sumstats = py_load_object(ref_table_sumstats_file, pickle = "pickle")
ref_table_latent_variables = py_load_object(ref_table_latent_variables_file, pickle = "pickle")

ref_table_sumstats = ref_table_sumstats[rownames(ref_table_latent_variables),]

columns_to_remove <- numeric()
for (i in seq_len(ncol(ref_table_sumstats))){
  if (sum(is.na(ref_table_sumstats[i])) > nrow(ref_table_sumstats)/10){
    columns_to_remove <- c(columns_to_remove,i) 
  }
}
ref_table_sumstats=ref_table_sumstats[-columns_to_remove]
rows_to_keep=complete.cases(ref_table_sumstats)
ref_table_sumstats=ref_table_sumstats[rows_to_keep,]
ref_table_latent_variables=ref_table_latent_variables[rows_to_keep,]

i=1

for (i in seq_len(ncol(ref_table_latent_variables))){
  param <- ref_table_latent_variables[i]
  reftable <- cbind(param,ref_table_sumstats)
  names(reftable) <- c("param",names(ref_table_sumstats))
  RFmodel <- regAbcrf(formula = log10(param)~.,
                      data = reftable,    
                      ntree = 1000,
                      paral = T)
  assign(paste0("rfmodel_",names(ref_table_latent_variables)[i]),RFmodel)
  
  
}


RFmodel <- regAbcrf(formula = log10(theta_10)~.,
                          data = reftable[c("theta_10",sumstats_names)],
                          ntree = 1000,
                          paral = T)
