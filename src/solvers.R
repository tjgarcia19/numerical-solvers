#***********************************************************************************************************************
#   CMSC 150 Numerical and Symbolic Computation Project
#   Created by Trestan Janos G. Garcia
#
#   This program contains the functions needed to run the application. It solves Polynomial Regression, Quadratic
#   Spline Interpolation and Simplex.
#
#   Date Created: December 6,2019
#***********************************************************************************************************************

#***********************************************************************************************************************
#This function reads the CSV file input and return the sorted x and y values in a list
loadCSV <- function(f){
  fileOpen= file(f, open = "rt")
  data = read.csv(fileOpen, header = FALSE)
  close(fileOpen)
  table = matrix(0, ncol = 2,nrow=length(data[,1]), dimnames = list(NULL,c("x","y")))
  table[,1]=data[,1];table[,2]=data[,2]
  #SORTING (It uses selection sort algo)
  i = 1;
  while(i<=length(data[,1])){
    j=i;
    while(j>1){
      if(table[j-1,1]>table[j,1]){
        temp = table[j-1,]      #swapping the rows
        table[j-1,] = table[j,]
        table[j,] = temp
      }else break;
      j=j-1
    }
    i=i+1
  }
  return(list(points = list(table[,1], table[,2]), xyTab = table))
}

#***********************************************************************************************************************
#This function computes the unknowns in a given set of equations, it accepts the number of equations
#and the augmented coefficient matrix, it returns the unknown values
GaussJordan2 <- function(numEq,matrix){
  pivotRow = 1
  #this loop is for moving the pivotrow from 1 to n
  while(pivotRow <= numEq){
    i = pivotRow; j = i; k = i;
    #checking the largest, before pivoting
    greatest = abs(matrix[pivotRow,j]) 
    while(i<=numEq){
      if(abs(matrix[i,j])>greatest){
        greatest = abs(matrix[i,j])
        k = i #this will determine the number of row that needed to be swapped to the current pivot row
      }
      i = i+1
    }
    #pivoting
    temp = c()
    temp = matrix[pivotRow, ]
    matrix[pivotRow, ] = matrix[k, ]
    matrix[k, ] = temp
    #NORMALIZE HERE
    if(matrix[pivotRow, pivotRow] != 0){
      matrix[pivotRow, ] = matrix[pivotRow, ] / matrix[pivotRow, pivotRow]
    }else{
      cat("\nNo Solution\n")
      return(NA)
    }
    #this loop is for elimination non-zero elements under/above the pivot element
    l = 1
    while(l<=numEq){
      if(l != pivotRow){ #if it is not the pivot row, eliminate non-zero
        pivotElement = matrix[pivotRow,j]
        valueToElim = matrix[l, j]
        multiplier = valueToElim/pivotElement
        vector1 = matrix[pivotRow, ] * multiplier
        vector2 = matrix[l, ] - vector1
        matrix[l, ] = vector2
      }
      l = l + 1
    }
    pivotRow = pivotRow + 1
  }
  #this loop adds the value of the unknowns in the solution vector
  result = c(); m=1;
  while(m<=numEq){
    result =c(result, matrix[m,numEq+1])
    m = m +1
  }
  return(result) #returning the solution vector
}

#***********************************************************************************************************************
#This function accepts an integer and a list. It returns a list with the augmented matrix, the values
#of unknowns, the polynomial string , and the polynomial function
PolynomialRegression <- function(degree, listOfPoints){
  if(degree <= 0){#this checks if the degree is greater than 0, if not, return NA
    cat("The degree must be greater than zero!\n")
    return(NA)
  }
  #getting the values from the list
  x = listOfPoints[[1]]
  y = listOfPoints[[2]]
  dimension = degree+1
  augMatrix = matrix(0, nrow = dimension, ncol= dimension+1)#creating an empty matrix
  
  i = 1;
  #this loop will put values in the matrix, and set up the augmented coefficient matrix
  while(i<=dimension){
    e = i - 1; eRHS = e ; j = 1;
    #e is the exponent
    #eRHS will store the first value of e to be used by the Right Hand Side
    while(j<=dimension+1){
      if(i == 1 && j == 1){ augMatrix[i,j] = length(x)
      }else if(j == dimension+1){ augMatrix[i,j] = sum((x^eRHS)*y)
      }else augMatrix[i,j] = sum(x^e);
      j = j + 1; e = e + 1;
    }
    i = i + 1
  }
  #this gets the values of unknowns using Gauss-Jordan Method
  unknowns = GaussJordan2(dimension, augMatrix)
  str = "function(x) " #initializing the string with the string function(x)
  k = 1; exp = 0 #exp will be the exponent
  #this loop concatenates the values of the unknowns to the string
  while(k<=dimension){
    if(k==1){#if it is the first no need to concat '*x^' since it will be raised to 0
      str = paste(str, toString(unknowns[k])," + ", sep = "")
    }else if(k == dimension){#if it is the last, dont add '+'
      str = paste(str, toString(unknowns[k]),"*x^",toString(exp),sep = "")
    }else{
      str = paste(str, toString(unknowns[k]),"*x^",toString(exp)," + ", sep = "")
    }
    k = k+1; exp = exp + 1;
  }
  #this makes the string a function
  str_function = eval(parse(text = str))
  
  return(list(augcoeffmatrix = augMatrix, unknowns = unknowns, polynomial_string = str, polynomial_function = str_function))
}

#***********************************************************************************************************************
#This function solves interpolation using Quadratic Splines. It accepts the number to be estimated and the list of
#x and y points. It returns a list containing the num of interval, functions per interval and the estimated value
QuadraticSpline <- function(num,listofPoints){
  x = listofPoints[[1]] #getting the x and y points
  y = listofPoints[[2]]
  numPoints = length(x)
  interval = numPoints-1
  augMatrix = matrix(0, nrow = 3*interval-1, ncol = 3*interval)#creating empty matrix
  table = matrix(0, nrow = numPoints, ncol = 2)
  table[,1] = x; table[,2] = y;#assigning the points of x and y to the 1st and 2nd col of the matrix
  
  if(length(x)==1 || num<table[1,1] || num>table[numPoints,1]){#this will only allow values in between the end points
    cat("The value cannot be estimated by the interpolating function.")
    return(NA)
  }
  
  #Condition 1: interior knots
  i=2; l=1;
  while(i<=interval){
    j=1
    while(j<=2){
      if(i==2 && j==1){ k=1; m=2;
      }else if(j==1){ k=2;
      }else if(j==2){ k=2; m=m+3 }
      while(k>=0){
        augMatrix[l,m-k] = table[i,1]^k;
        k=k-1;
      }
      augMatrix[l,3*interval]=table[i,2]
      j=j+1;l=l+1
    }
    i=i+1
  }
  
  #Condition 2: end knots
  augMatrix[l,1]=table[1,1]
  augMatrix[l,2]=1
  augMatrix[l, 3*interval]=table[1,2]
  l=l+1
  i=3
  while(i>0){
    augMatrix[l,3*interval-i]=table[numPoints,1]^(i-1)
    i=i-1
  }
  augMatrix[l,3*interval]=table[numPoints,2]
  l=l+1
  
  #Condition 3: first derivative of the interior knots
  i=2;j=1
  while(i<=interval){
    if(i==2){
      augMatrix[l,j]=1
      augMatrix[l,j+2]=-table[i,1]*2
      augMatrix[l,j+3]=-1
      j=j+2
    }else{
      augMatrix[l,j]=table[i,1]*2
      augMatrix[l,j+1]=1
      augMatrix[l,j+3]=-table[i,1]*2
      augMatrix[l,j+4]=-1
      j=j+3
    }
    i=i+1;l=l+1
  }
  
  unknowns = GaussJordan2(3*interval-1,augMatrix) #solving the unknowns
  
  #creating the list of function per interval
  i=1
  numVar = length(unknowns)
  tempEqs2 = c()
  str2 = "function(x) "
  k=1
  while(i<=interval){
    if(i==1){
      eq = paste(str2, toString(unknowns[k]),"*x+", toString(unknowns[k+1]), sep="")
      k=k+2
    }else{
      eq = paste(str2, toString(unknowns[k]),"*x^2+", toString(unknowns[k+1]), "*x+", toString(unknowns[k+2]), sep="")
      k=k+3
    }
    i=i+1
    tempEqs2 = c(tempEqs2, eq)
  }
  #this will determine which function in the interval will be used to estimate the value
  if(table[numPoints,1]==num){ i=interval;
  }else if(table[numPoints,1]==num){ i=1;
  }else{
    i=numPoints;
    while(table[i,1]>num){i=i-1;}
  }
  
  est_fxn = eval(parse(text=tempEqs2[i]))#makes the string function
  est = est_fxn(num)

  return(list(num_interval=interval,interval_fxns = tempEqs2, estimating_fxn=est_fxn, estimate=est))
}

#***********************************************************************************************************************
#This function finds the minimum cost to ship materials from one plant to another. It accepts a vector of 
#cost, supply and demand. It returns a list containing the initial tableau and per iteration, the solution table and the number of iteration
Simplex <- function(cost, supply, demand){
  augMatrix = matrix(0, nrow = 16, ncol=25, dimnames = list(c(1:16), c("y1", "y2","y3","y4","y5","y6","y7","y8","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","Z","Sol'n"))) #
  
  #this sets up the initial tableau
  j=1; i=1; k=1; m=1; demandCnt=1; supplyCnt=1;
  while(j<=25){
    if(j<=5){
      augMatrix[i,j]=1
      augMatrix[i+5,j]=1
      augMatrix[i+10,j]=1
      augMatrix[16,j]=demand[demandCnt]*-1
      i=i+1;demandCnt=demandCnt+1;
    }else if(j<=8 && j>5){
      l=1
      while(l<=5){
        augMatrix[k,j]=-1
        l=l+1;k=k+1
      }
      augMatrix[16,j] = supply[supplyCnt]
      supplyCnt=supplyCnt+1
    }else if(j<=24 && j>8){
      augMatrix[m,j]=1
      m=m+1
    }else if(j==25){
      i=1
      while(i<=15){
        augMatrix[i,j]=cost[i]
        i=i+1
      }
    }
    j=j+1
  }
  init_tableau = augMatrix
  
  #Compute here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iteration=0; isOptimized = FALSE
  tableauPerIteration=list(); solnPerIteration=list();
  while(iteration!=1000){#max iteration is 1000 if it exceeds, the problem is assumed to have no feasible solution
    count=0
    smallest = augMatrix[16,1]
    j=1; pivotCol = j;
    while(j<25){#this loop finds the pivot column
      if(augMatrix[16,j] < smallest){
        smallest=augMatrix[16,j]
        pivotCol = j
      }
      if(augMatrix[16,j]>=0) count=count+1;#this counts the negative in the bottom row
      j=j+1
    }

    if(count==24){isOptimized=TRUE}#if count is 24, it is optimized
    if(isOptimized){break}
    
    testRatio=999999
    i=1; pivotRow = i
    while(i<16){#this finds the pivot row
      if(augMatrix[i,pivotCol] > 0){
        sampRatio = augMatrix[i,25]/augMatrix[i,pivotCol]
        if(sampRatio < testRatio){
          testRatio=sampRatio;
          pivotRow = i
        }
      }
      i=i+1
    }
    #NORMALIZE HERE, if the denominator is 0, stop the iteration, no feasible solution
    if(augMatrix[pivotRow, pivotCol]!=0){ augMatrix[pivotRow, ] = augMatrix[pivotRow, ] / augMatrix[pivotRow, pivotCol]
    }else{break}
    
    #this loop is for elimination non-zero elements under/above the pivot element
    l = 1
    while(l<=16){
      if(l != pivotRow && augMatrix[l,pivotCol]!=0){ #if it is not the pivot row, eliminate non-zero
        pivotElement = augMatrix[pivotRow,pivotCol]
        valueToElim = augMatrix[l, pivotCol]
        multiplier = valueToElim/pivotElement
        vector1 = augMatrix[pivotRow, ] * multiplier
        vector2 = augMatrix[l, ] - vector1
        augMatrix[l, ] = vector2
      }
      l = l + 1
    }
    
    tableauPerIteration[[toString(iteration)]]=augMatrix
    
    #Setup answer here
    unknowns = augMatrix[16,9:23]
    solutionTable = matrix(0, nrow=4, ncol=6, dimnames = list(c("Denver", "Phoenix", "Dallas", "Total"), c("Total", "Sacramento", "Salt Lake", "Albuquerque", "Chicago", "New York")))
    i=1;k=1
    while(i<=3){
      j=2
      while(j<=6){
        solutionTable[i,j]=unknowns[k]
        k=k+1;j=j+1
      }
      i=i+1
    }
    j=2
    while(j<=6){
      solutionTable[4,j]=sum(solutionTable[1:3,j])
      j=j+1
    }
    i=1;k=1
    while(i<=3){
      j=2; total = 0
      while(j<=6){
        total = total + solutionTable[i, j]*cost[k]
        j=j+1;k=k+1
      }
      solutionTable[i,1]=total
      i=i+1;j=j+5
    }
    solutionTable[4,1]=sum(solutionTable[1:3,1])
    solnPerIteration[[toString(iteration)]]=solutionTable
    iteration=iteration+1
  }

  return(list(init_tableau=init_tableau,num_iteration = iteration,tableau_iteration=tableauPerIteration,soln_iteration=solnPerIteration,solution=solutionTable,optimize=isOptimized))
}
