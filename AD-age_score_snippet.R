# Snippet to compute the AD-age scores 
# considering two types of weight(age) functions

control.score <- function(weight.age){
  return (log(1-weight.age)-0.5)
}

case.score <- function(weight.age){
  return (-log(weight.age)+0.5)
}

# linear age weight
linear.weight<-function(age){
  age <- as.numeric(age)
  w = (age-59.5)/(100.5-59.5)
  return (w)
}

# piecewise age weight
piecewise.weight<-function(age){
  age <- as.numeric(age)
  RANGE = 32 
  if (age < 60) {
    w <- 5/(10*RANGE)
  } else if (age < 65) {
    w <- (age-55)/(10*RANGE)
  } else if (age < 75) {
    w <- 4*(age-55)/(10*RANGE) - 3/RANGE
  } else if (age < 80) {
    w <- 10*(age-55)/(10*RANGE) - 15/RANGE
  } else if (age < 90) {
    w <- 16*(age-55)/(10*RANGE) - 30/RANGE
  } else if (age < 99) {
    w <- 6*(age-55)/(10*RANGE) + 5/RANGE
  } else if (age >= 99) {
    w <- 0.975
  }
  return (w)
}

piecewise.score <- function(row){
  if (row['diagnosis'] == 'CN'){
    return (control.score(piecewise.weight(row['age'])))
  }else if (row['diagnosis'] == 'AD'){
    return (case.score(piecewise.weight(row['age'])))
  }else{
    return (NA)
  }
}

linear.score <- function(row){
  if (row['diagnosis'] == 'CN'){
    return (control.score(linear.weight(row['age'])))
  }else if (row['diagnosis'] == 'AD'){
    return (case.score(linear.weight(row['age'])))
  }else{
    return (NA)
  }
}

# age range
age.range <- seq(60, 100, by=0.5)

# here you can load your own dataframe instead of generating random data
age <- sample(age.range, replace=T, size=1000)
diagnosis <- sample(c('CN', 'AD'), replace=T, size=1000)
df <- as.data.frame(cbind(age, diagnosis))
df$age <- as.numeric(df$age)

df[,'linear.score'] <- apply(df, 1, linear.score)
df[,'piecewise.score'] <- apply(df, 1, piecewise.score)
