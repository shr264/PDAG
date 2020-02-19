coords4 <- function(){
  coords = matrix(0,ncol = 2, nrow = 30)
  num2 = 0
  for (i in 1:6){
    coords[i,2] = ifelse(i%%2==0,7,8)
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = 0
  for (i in 7:16){
    coords[i,2] = ifelse(i%%2==0,5,6)
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0
  for (i in 17:23){
    coords[i,2] = ifelse(i%%2==0,3,4)
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0
  for (i in 24:30){
    coords[i,2] = ifelse(i%%2==0,1,2)
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords5 <- function(){
  coords = matrix(0,ncol = 2, nrow = 30)
  num2 = 0
  for (i in 1:6){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 7:14){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 15:21){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 22:28){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 29:30){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords5wTI <- function(){
  coords = matrix(0,ncol = 2, nrow = 25)
  num2 = 0
  for (i in 1:6){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 7:11){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 12:18){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 19:23){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 24:25){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords5woTI <- function(){
  coords = matrix(0,ncol = 2, nrow = 24)
  num2 = 0
  for (i in 1:6){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 7:11){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 12:18){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 19:22){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 23:24){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords8 <- function(){
  coords = matrix(0,ncol = 2, nrow = 30)
  num2 = 0
  for (i in 1:2){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 3:3){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 4:4){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 11:14){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 15:21){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.5
  for (i in 22:23){
    coords[i,2] = 3
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.5
  for (i in 24:28){
    coords[i,2] = 2
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0
  for (i in 29:30){
    coords[i,2] = 1
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords10 <- function(){
  coords = matrix(0,ncol = 2, nrow = 30)
  num2 = 0
  for (i in 1:2){
    coords[i,2] = 10
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 3:3){
    coords[i,2] = 9
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 4:4){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 5:6){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 7:10){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 11:14){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 15:21){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.5
  for (i in 22:23){
    coords[i,2] = 3
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.5
  for (i in 24:28){
    coords[i,2] = 2
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0
  for (i in 29:30){
    coords[i,2] = 1
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords9wTI <- function(){
  coords = matrix(0,ncol = 2, nrow = 25)
  num2 = 0
  for (i in 1:2){
    coords[i,2] = 10
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 3:3){
    coords[i,2] = 9
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 4:4){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 5:6){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 7:9){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 10:11){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 12:18){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.5
  for (i in 19:23){
    coords[i,2] = 3
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.5
  for (i in 24:25){
    coords[i,2] = 2
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

coords9woTI <- function(){
  coords = matrix(0,ncol = 2, nrow = 24)
  num2 = 0
  for (i in 1:2){
    coords[i,2] = 10
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  num2 = -0.75
  for (i in 3:3){
    coords[i,2] = 9
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 4:4){
    coords[i,2] = 8
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.75
  for (i in 5:6){
    coords[i,2] = 7
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.75
  for (i in 7:9){
    coords[i,2] = 6
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -1
  for (i in 10:11){
    coords[i,2] = 5
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 1
  for (i in 12:18){
    coords[i,2] = 4
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = -0.5
  for (i in 19:22){
    coords[i,2] = 3
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  
  num2 = 0.5
  for (i in 23:24){
    coords[i,2] = 2
    coords[i,1] = 2*num2
    num2 = num2 + 1
  }
  return(coords)
}

get_coords <- function(num_groups){
  if(num_groups==8){coords = coords8()
  } else if(num_groups==4){coords = coords4()
  } else if(num_groups==5){coords = coords5()
  } else if(num_groups==10){coords = coords10()}
  return(coords)
}