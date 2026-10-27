## Aims to summary the min distance from cell type1 to cell type2 (or spots set) in spatial transcriptome.

## The input need the two cell-barcode sets: barcode_set.1,barcode_set.2, which indicate the interseted cell types.
## Also need to provide the coordinates of all spots.
## the cell type information is stored in object@meta.data

## example
# tmp_object <- subset(KRAS_Pt8_Tm_result,subset = layer == 'Mucosal_layer')
# coord <- tmp_object@images$sliceD1@coordinates
# weight.matr <- dplyr::select(tmp_object@meta.data,'CTHRC1.Stromal5',Mal)
# dist <- summ_dist_sampling(object = tmp_object,
#                      barcode_set.2 = colnames(tmp_object)[which(tmp_object$subcluster %in% c(9,3,9.2))],
#                      barcode_set.1 = colnames(tmp_object)[which(tmp_object@meta.data[['CTHRC1.Stromal5']] > 0 )],
#                      coord,weight = T,weight.matr = weight.matr,sampling.num = 1)
# dist <- summ_dist_V2(object = tmp_object,coord = coord,
#                            barcode_set.2 = colnames(tmp_object)[which(tmp_object$subcluster %in% c(9,3,9.2))],
#                            barcode_set.1 = colnames(tmp_object)[which(tmp_object@meta.data[['CTHRC1.Stromal5']] > 0 )],
#                            weight = T,weight.matr = weight.matr)
# dist
# 
# 
# # 2
# table(KRAS_Pt8_Tm_result$layer)
# layer_list <- c('Muscular_layer','Mucosal_layer','Epithelial_layer') %>% as.list()
# stromal_name_list <- paste(names(mk_score_list_overlap),1:5,sep = '') %>% as.list()
# 
# dist_stromal_layer_list <- lapply(layer_list, function(x){
#   stromal_name_list <- paste(names(mk_score_list_overlap),1:5,sep = '') %>% as.list()
#   tmp_object <- subset(KRAS_Pt8_Tm_result,subset = layer == x)
#   coord <- tmp_object@images$sliceD1@coordinates
#   dist_stromal <- lapply(stromal_name_list, function(y){
#     if(length(colnames(tmp_object)[which(tmp_object@meta.data[[y]] > 0 )]) < 0.1*ncol(tmp_object))
#       return(-2)
#     weight.matr <- dplyr::select(tmp_object@meta.data,y,Mal)
#     dist <- summ_dist_V2(object = tmp_object,
#                          barcode_set.2 = colnames(tmp_object)[which(tmp_object$subcluster %in% c(9,3,9.2))],
#                          barcode_set.1 = colnames(tmp_object)[which(tmp_object@meta.data[[y]] > 0 )],
#                          coord,weight = 'double',weight.matr = weight.matr)
#     dist
#   })
#   dist_stromal
# })
# dist_stromal_layer_list
# 
# 
# layer_list <- layer_list[[2]]
# tmp_object
# 
# avg_stromal_layer_list <- lapply(layer_list, function(x){
#   stromal_name_list <- paste(names(mk_score_list_overlap),1:5,sep = '') %>% as.list()
#   tmp_object <- subset(KRAS_Pt8_Tm_result,subset = layer == x)
#   coord <- tmp_object@images$sliceD1@coordinates
#   dist_stromal <- lapply(stromal_name_list, function(y){
#     if(length(colnames(tmp_object)[which(tmp_object@meta.data[[y]] > 0 )]) < 0.1*ncol(tmp_object))
#       return(-2)
#     weight.matr <- dplyr::select(tmp_object@meta.data,y,Mal)
#     dist <- summ_dist_sampling(object = tmp_object,
#                          barcode_set.2 = colnames(tmp_object)[which(tmp_object$subcluster %in% c(9,3,9.2))],
#                          barcode_set.1 = colnames(tmp_object)[which(tmp_object@meta.data[[y]] > 0 )],
#                          coord,weight = 'double',weight.matr = weight.matr,sampling.num = 10,sampling.size = 0.1*ncol(tmp_object))
#     dist
#   })
#   dist_stromal
# })
# avg_stromal_layer_list


# dist_stromal_layer_list <- lapply(layer_list, function(x){
#   stromal_name_list <- paste(names(mk_score_list_overlap),1:5,sep = '') %>% as.list()
#   tmp_object <- subset(KRAS_Pt8_Tm_result,subset = layer == x)
#   coord <- tmp_object@images$sliceD1@coordinates
#   dist_stromal <- lapply(stromal_name_list, function(y){
#     weight.matr <- dplyr::select(tmp_object@meta.data,y,Mal)
#     dist <- summ_dist_V2(object = tmp_object,
#                          barcode_set.2 = colnames(tmp_object)[which(tmp_object$subcluster %in% c(9,3,9.2))],
#                          barcode_set.1 = colnames(tmp_object)[which(tmp_object@meta.data[[y]] > 0 )],
#                          coord,weight = 'double',weight.matr = weight.matr)
#     dist
#   })
#   dist_stromal
# })

summ_dist_V2 <- function(object,barcode_set.1,barcode_set.2,coord,weight = 'double',weight.matr){
  coord$tag1 <- 0
  coord$tag2 <- 0
  coord$tag1[which(rownames(coord) %in% barcode_set.1)] <- 1
  coord$tag2[which(rownames(coord) %in% barcode_set.2)] <- 2
  # if(weight != F)
  #   if(is.na(weight.matr[1,1]))
  #     stop("'weight.matr' needs to be provided if you want to take weights into account.")
  #   else 
  #     weight.matr <- matrix(data = 1,nrow = nrow(coord),ncol = 2)
  
  # Use BFS to calculate the min dist
  BFS_search <- function(x,y)
  {
    coord$tag <- 0
    index_a <- array(dim = 2000)
    dist_a = array(dim = 2000)
    index_a[1] <- which(coord$row == x & coord$col == y)
    dist_a[1] = 0
    head = 1
    tail = 2
    
    existornot <- function(x,y){
      if(sum(which(coord$row == x & coord$col == y)) == 0)
        return(0)
      else
        if(coord$tag[which(coord$row == x & coord$col == y)] == 1)
          return(0)
      return(1)
    }
    
    while(head < tail) {
      if (coord$tag2[index_a[head]] == 2)
      {
        if(weight == F)
          return(dist_a[head])
        else 
          if(weight == 'single')
            return(list(dist_a[head],weight.matr[index_a[1]]))
        else 
          if(weight == 'double')
          {
            dist_tmp <- dist_a[head]
            weight.max <- weight.matr[index_a[head],2]
            while(dist_a[head]== dist_tmp & head < tail){
              weight.max <- max(weight.max,weight.matr[index_a[head],2])
              head <- head + 1
            }
            return(list(dist_tmp,weight.matr[index_a[1],1]*weight.max))
          }
        
      }
      coord$tag[index_a[head]] <- 1
      
      if(existornot(coord$row[index_a[head]]+1,coord$col[index_a[head]]-1))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]]+1 & coord$col == coord$col[index_a[head]]-1)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      if(existornot(coord$row[index_a[head]]-1,coord$col[index_a[head]]-1))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]]-1 & coord$col == coord$col[index_a[head]]-1)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      if(existornot(coord$row[index_a[head]],coord$col[index_a[head]]-2))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]] & coord$col == coord$col[index_a[head]]-2)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      if(existornot(coord$row[index_a[head]],coord$col[index_a[head]]+2))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]] & coord$col == coord$col[index_a[head]]+2)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      if(existornot(coord$row[index_a[head]]+1,coord$col[index_a[head]]+1))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]]+1 & coord$col == coord$col[index_a[head]]+1)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      if(existornot(coord$row[index_a[head]]-1,coord$col[index_a[head]]+1))
      {
        index_a[tail] <- which(coord$row == coord$row[index_a[head]]-1 & coord$col == coord$col[index_a[head]]+1)
        coord$tag[index_a[tail]] <- 1
        dist_a[tail] <- dist_a[head]+1
        tail <- tail+1
      }
      
      
      head <- head + 1
    }
    return(list(-1,
                weight.matr[index_a[1],1]*mean_dist_2))
  }
  
  min_dist_a <- array(dim = length(barcode_set.1))
  weight_a <- array(dim = length(barcode_set.1))
  tmp_list <- list(x=0,y=0)
  coord_tag1 <- dplyr::filter(coord,tag1 == 1)
  mean_dist_2 <- mean(weight.matr[which(rownames(weight.matr) %in% barcode_set.2),2])
  pb <- txtProgressBar(style=3)
  for (i in 1:length(barcode_set.1)) {
    tmp_list <- BFS_search(coord_tag1$row[i],coord_tag1$col[i])
    min_dist_a[i] <- tmp_list[[1]]
    weight_a[i] <- tmp_list[[2]]
    setTxtProgressBar(pb, i/length(barcode_set.1))
  }
  #summ_result <- cbind(unlist(min_dist_a),unlist(weight_a)) %>% as.data.frame()
  #colnames(summ_result) <- c("type",'weight')
  #summ_result <- dplyr::group_by(summ_result,type) %>%  dplyr::summarise(sum = sum(weight))
  close(pb)
  return(list(min_dist_a,weight_a))
}


summ_dist_sampling <- function(object,barcode_set.1,barcode_set.2,coord,sampling.size = 100, sampling.num = 1000,weight = T,weight.matr = NULL){
  if(weight != F)
  {
    if(is.na(weight.matr))
      stop("'weight.matr' needs to be provided if you want to take weights into account.")
  }
  else 
    weight.matr <- matrix(data = 1,nrow = nrow(coord),ncol = 2)
  
  time_start <- proc.time()
  # start sampling
  sampling_a <- array(dim = sampling.num)
  for (i in 1:sampling.num)
  {
    # randomly select several spots from barcode_set.1, size = sampling.size
    set.seed(i)
    spot_select <- sample(barcode_set.1,size = min(length(barcode_set.1),sampling.size),)
    # run summ_dist
    dist <- summ_dist_V2(object,barcode_set.1 =spot_select, barcode_set.2 = barcode_set.2,coord = coord,weight = weight,weight.matr = weight.matr)
    # record the avg_dis
    dist$percentage.adj <- dist$sum / sum(dist$sum)
    avg_value <- dplyr::filter(dist,type != -1) %>% 
      dplyr::summarise(.,value = sum(percentage.adj*type)) 
    sampling_a[i] <- avg_value
  }
  time_end <- proc.time()
  running_time <- time_end - time_start
  cat(running_time[[3]][1])
  return(sampling_a)
}

# create permutation test
permutation_test <- function(object,n.spot,barcode_set.2,coord,test.num = 1000,weight = F,weight.matr = NULL){
  permu_a <- array(dim = test.num)
  if(weight == 'single')
    colnames(weight.matr) <- 'set.1'
  if(weight == 'double')
    colnames(weight.matr) <- c('set.1','set.2')
  for (i in 1:test.num) {
    set.seed(i)
    spot_select <- sample(colnames(object),size = n.spot)
    if(weight!= F)
    {
      weight.new <- dplyr::filter(weight.matr,set.1 > 0)
      weight.matr$set.1 <- 0
      weight.matr[spot_select,1] <- weight.new
    }
    
    dist <- summ_dist_V2(object,barcode_set.1 =spot_select, barcode_set.2 = barcode_set.2,coord = coord,weight = weight,weight.matr = weight.matr)
    dist$percentage.adj <- dist$sum / sum(dist$sum)
    avg_value <- dplyr::filter(dist,type != -1) %>% 
      dplyr::summarise(.,value = sum(percentage.adj*type)) 
    permu_a[i] <- avg_value
  }
  return(permu_a)
}
