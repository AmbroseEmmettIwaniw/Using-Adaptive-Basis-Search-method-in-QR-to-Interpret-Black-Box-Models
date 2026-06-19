normalize_data <- function(data, response_index){
  Response <- data[,response_index]
  data_norm <- data[,-response_index]
  
  min_data <- as.double(apply(data_norm,2,min))
  max_data <- as.double(apply(data_norm,2,max))
  
  data_norm = sweep(data_norm, 2, min_data)
  data_norm = sweep(data_norm, 2, max_data - min_data, FUN='/')
  
  data_norm = cbind(data_norm, Response)
  
  data_norm
}