check_create_dir <- function(new_dir){
  # check if directory exists. If it doesn't create it. 
  if (!dir.exists(new_dir)){
    dir.create(new_dir)
  }
}