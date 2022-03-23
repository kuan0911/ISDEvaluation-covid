aveMetrics = function(ISD) {
  if(length(ISD$survivalCurves)>2) {print('aveMetrics can not work with more than one model')}
  print(paste(ISD$results$Model[1],' N =',ISD$results$N[1],sep=' '))
  C = mean(ISD$results$Concordance)
  Cstd = sd(ISD$results$Concordance)
  print(paste('Conc:',trunc(C*1000)/1000,'+/-',trunc(Cstd*1000)/1000,sep=' '))
  B = mean(ISD$results$BrierInt)
  Bstd = sd(ISD$results$BrierInt)
  print(paste('Brier:',trunc(B*1000)/1000,'+/-',trunc(Bstd*1000)/1000,sep=' '))
  L = mean(ISD$results$L1Loss)
  Lstd = sd(ISD$results$L1Loss)
  print(paste('L1Loss:',trunc(L*1000)/1000,'+/-',trunc(Lstd*1000)/1000,sep=' '))
  D = ISD$results$DCalibration[1]
  print(paste('DCal:',trunc(D*1000)/1000,sep=' '))
  return(list(Concordance=C,Concordance_std=Cstd,BrierInt=B,BrierInt_std=Bstd,L1Loss=L,L1Loss_std=Lstd,DCalibration=D))
}