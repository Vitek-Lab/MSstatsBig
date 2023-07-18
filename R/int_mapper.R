#' @export
string_to_int_mapper = function(data, string_cols, int_cols){

  data = as.data.table(data)
  id_cols = paste0("id_", string_cols)

  data[, (string_cols) := lapply(.SD, as.factor), .SDcols=string_cols]
  data[, (id_cols) := lapply(.SD, as.numeric), .SDcols=string_cols]

  lookup_list = list()
  for (i in seq_along(string_cols)){
    lookup_cols = c(string_cols[[i]], id_cols[[i]])
    lookup_list[[string_cols[[i]]]] = unique(data[, ..lookup_cols])
  }

  cols = c(id_cols, int_cols)
  data = data[,..cols]

  setnames(data, id_cols, string_cols)

  return(list(data=data, lookup=lookup_list))
}

#' @export
int_to_string_mapper = function(data, lookup_list, string_cols){

  for (i in seq_along(lookup_list)){
    data = merge(data, lookup_list[[i]], all.x=TRUE,
                 by.x = names(lookup_list)[[i]],
                 by.y=paste0("id_", names(lookup_list)[[i]]))
  }
  data[, (string_cols):=NULL]

  setnames(data, paste0(string_cols, ".y"), string_cols)
  return(data)
}
