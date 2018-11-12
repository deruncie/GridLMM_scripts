library(GridLMM)
library(pryr)
library(foreach)
library(cowplot)
library(reshape2)

ns = floor(2^seq(8,12,by=.5))

times = foreach(n = ns,.combine = 'rbind') %do% {
  print(n)
  K = tcrossprod(matrix(rnorm(n^2),n)) + diag(1,n)
  X = matrix(rnorm(n*1e4),n)
  
  t_chol = as.numeric(system.time(chol(K))[3])
  t_KX = as.numeric(system.time(K %*% X)[3])
  
  s_K = as.numeric(object_size(K))
  c(n,t_chol,t_KX,s_K)
}
write.csv(times,file = 'Results/Sample_size/sample_size_times.csv',row.names=F)

time_results = read.csv('Results/Sample_size/sample_size_times.csv')
colnames(time_results) = c('N','Time_cholesky','Time_RX','Bytes_cholesky')

time_results$Time_RX1 = time_results$Time_RX / 1e4   # per marker, only upper-triangle
time_results$Bytes_cholesky = time_results$Bytes_cholesky # only store upper-triangle

model_Time_cholesky = lm(Time_cholesky^(1/3) ~ N,time_results) # model as linear in n^3
model_Time_RX = lm(Time_RX1^(1/2) ~ N,time_results) # model as linear in n^2


ns = floor(2^seq(8,log2(1e5),length=20))
fitted_times = data.frame(N = ns,
                          `Cholesky (nxn)` = predict(model_Time_cholesky,list(N=ns))^3,
                          RX1 = predict(model_Time_RX,list(N=ns))^2,check.names = FALSE)
fitted_times$`RX (nx1e5)` = 1e5*fitted_times$RX1
fitted_times$Total_GridLMM = fitted_times$Cholesky + fitted_times$`RX (nx1e5)`
fitted_times$Total_LDAK = 1e5 * (fitted_times$Cholesky + fitted_times$`RX (nx1e5)`)

fitted_times_tall = melt(fitted_times,id.vars = 'N')
fitted_times_tall = subset(fitted_times_tall,variable != 'RX1')

p1 = ggplot(fitted_times_tall,aes(x=N)) + 
  geom_line(aes(y=value,color = variable),size=1) + 
  geom_point(data = time_results,aes(y=Time_cholesky),col='red',size=2) + 
  geom_point(data = time_results,aes(y=1e5*Time_RX1),col='darkgreen',size=2) + 
  scale_x_log10() + scale_y_log10() + 
  background_grid(major = "xy", minor = "xy") +
  theme(legend.title = element_blank()) + 
  ylab('seconds per iteration') + xlab('N observations')
p1

model_Bytes_cholesky = lm(Bytes_cholesky^(1/2)~N,time_results)
fitted_size = data.frame(N = ns,
                          `MB_Cholesky` = predict(model_Bytes_cholesky,list(N=ns))^2/2e6
                          )
p2 = ggplot(fitted_size,aes(x=N)) + 
  geom_line(aes(y=MB_Cholesky),size=1) + 
  geom_point(data = time_results,aes(y=Bytes_cholesky/2e6),size=2)+
  scale_x_log10() + scale_y_log10() + 
  background_grid(major = "xy", minor = "xy") +
  theme(legend.title = element_blank()) + 
  ylab('MB for Cholesky matrix') + xlab('N observations')
p2

ss_time_size = plot_grid(p2,p1,nrow = 1,labels = c('a','b'),rel_widths = c(1,2));ss_time_size
save_plot('Figures/Sample_size_times.pdf',ss_time_size,base_aspect_ratio = 2.5)

