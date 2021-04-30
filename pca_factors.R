library(factoextra)
p=15
x1=c(29.60,25.00,24.33,24.02,22.00,17.64,11.82,6.90,6.70) #poverty rate
x2=c(69.3, 70.29, 71, 71.43, 71.7, 72, 72.72, 74.37, 77.77) #literacy rate
x3=c(68.72, 68.37, 68, 67.62, 67.22, 66.82, 66.4, 65.97, 65.53) #rural population%
x4=c(38.03, 38.51, 38.96, 39.41, 39.86, 40.29, 40.72, 41.15, 41.57) #population density
x5=c(5.64, 5.65, 5.67, 5.61, 5.57, 5.51, 5.42, 5.33, 5.36) #unemployment rate
x6=c(63.46, 65.54, 68.57, 72.81, 77.66, 83, 87.33, 92.09, 94.95) #per capita income
x7=c(4.18, 6.12, 7.04, 8.26, 8.00, 7.41, 6.39, 5.46, 5.24) #economic growth rate
x8=c(3.54, 3.54, 3.51, 3.60, 3.62, 3.75, 3.33, 3.25, 3.27) #healthcare spending(%of gdp)
x9=c(7.66, 4.86, 2.49, 4.94, 5.87, 6.35, 10.91, 9.31, 8.86) #inflation rate
x10=c(27.00, 27.00, 28.60, 28.80, 30.6, 30.60, 32.60, 32.60, 34.80) #percentage of smokers
x11=c(27.39, 27.64, 28.48, 29.36, 30.27, 31.22, 32.20, 33.22, 34.77) #participation rate (economically active)
x12=c(14.00, 14.20, 14.40, 14.70, 15.30, 15.90, 16.30, 16.30,16.50) #hunger statistics
x13=c(100000000, 137000000, 190099500, 243199000, 350000000, 375000000, 462000000, 462000000, 560000000) #internet users
x14=c(38045000, 60545120, 83272560, 106000000, 134000000, 136000000, 191000000, 250000000, 310000000) #social media users
x15=c(858370000, 934100000, 910200000, 886300000, 896300000, 1012000000, 1059000000, 1185360000, 1190000000) #mobile users

#scaling of data
x=matrix(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15),nrow=9,ncol=15,byrow=FALSE)     #data matrix
x=abs(scale(x))

#covariance matrix
co=cov(x)

#testing: lawley procedure
cm=cor(x)
i=1
ri=c(1:p)*0
for(i in 1:p)
{
ri[i]=c((sum(cm[,i])-1)/(p-1))
}
rj=(sum(cm)-15)/(p*(p-1)) #r bar

s2=0
for(i in 1:p)
{
for(j in 1:p)
{
if(i!=j)
{
s2=s2+((cm[i,j]-rj)*(cm[i,j]-rj))
}
}
}
s3=s2/2   #summation 1

s4=0
for(i in 1:15)
{
s4=s4+((ri[i]-rj)^2)   #summation 2
}

ga=(((p-1)^2)*(1-(1-rj)^2))/(p-(p-2)*((1-rj)^2))  #gamma

n=9
T=((n-1)/((1-rj)^2))*(s3-(ga*s4))
T
T_tab=128.804 #chi-square value at 5% los with ((p+1)(p-2))/2=104 dof 
if(T>T_tab){
print(' we reject our null hypothesis, hence all the eigen values are distinct and pca is applicable')
}else{
print('pca cannot be applied')
}

#eigenvalues and eigenvectors
ei=eigen(co)
eva=ei$values
eve=ei$vectors

#principal components
result=t(eve)%*%t(x)

#variation explained by principal components
pc.pca=prcomp(x, scale=FALSE)
summary(pc.pca)
fviz_eig(pc.pca)

#principal components to be chosen
print('Principal Components that effectively summarises total variance are:')
for(i in 1:p)
{
if(eva[i]>mean(eva))
{
print(i)
w=i
}
}
print('Their value is represented by each column of following vector:')
result[,1:w]

#correlation between original variables x and principal components y
cpc=matrix(nrow=15,ncol=9,byrow=TRUE,dimnames=list(c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15"),c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9")))
for(i in 1:p)
{
for(j in 1:n)
{
cpc[i,j]=((sqrt(eva[j]))*(eve[i,j]))/(sqrt(co[i,i]))
}
}
print(cpc)
