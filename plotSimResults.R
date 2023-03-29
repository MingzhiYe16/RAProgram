
#png(filename = '~/Desktop/fig2.png', width = 1000, height = 500)

boxplot(results$DistanceListFromEM ~ results$sampSize,
        at = c(1,4,7,10,13,16,19,22,25), border = brewer.pal(9,'Set1')[5],
        ylab='Total Variation Distance from Truth', xlab='Cells', main='Simulation Accuracy',
        axes = FALSE, ylim=c(0,max(c(results$DistanceListFromEM,results$DistanceListFromUN))), pch=16)
boxplot(results$DistanceListFromUN ~ results$sampSize,
        at = c(2,5,8,11,14,17,20,23,26), border = brewer.pal(9,'Set1')[4],
        axes = FALSE, add = TRUE, pch=16)

axis(2)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5), labels=c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4))
abline(v=c(3,6,9,12,15,18,21,24), lty=2, col='grey50')
legend('topright', fill=brewer.pal(9,'Set1')[4:5], border=NA, legend=c('Unique','EM'))






boxplot(results$KLDivergenceListFromEM ~ results$sampSize,
        at = c(1,4,7,10,13,16,19,22,25), border = brewer.pal(9,'Set1')[5],
        ylab='Kullback-Leibler Divergence from Truth', xlab='Cells', main='Simulation Accuracy',
        axes = FALSE, ylim=c(0,max(c(results$KLDivergenceListFromEM,results$KLDivergenceListFromUN))), pch=16)
boxplot(results$KLDivergenceListFromUN ~ results$sampSize,
        at = c(2,5,8,11,14,17,20,23,26), border = brewer.pal(9,'Set1')[4],
        axes = FALSE, add = TRUE, pch=16)

axis(2)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5), labels=c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4))
abline(v=c(3,6,9,12,15,18,21,24), lty=2, col='grey50')
legend('topright', fill=brewer.pal(9,'Set1')[4:5], border=NA, legend=c('Unique','EM'))




# for comparing simulation to real data

boxplot(results$DistanceListUN_EM ~ results$sampSize,
        at = c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4), log = 'x',
        xlim = c(70,7e4), boxwex = .6,
        xlab = '# Cells', ylab = 'Distance between EM and Unique')











###########
### OLD ###
###########


medians <- cbind(rowMedians(d100),rowMedians(d200),rowMedians(d500),rowMedians(d1K),rowMedians(d2K),rowMedians(d5K),rowMedians(d10K),rowMedians(d20K),rowMedians(d50K))


# total variation distance
###########################

png(filename = '~/Desktop/fig2.png', width = 1000, height = 500)

layout(matrix(1:2,nrow=1))

boxplot(cbind(t(d100)[,1:2],t(d200)[,1:2],t(d500)[,1:2],t(d1K)[,1:2],t(d2K)[,1:2],t(d5K)[,1:2],t(d10K)[,1:2],t(d20K)[,1:2],t(d50K)[,1:2]),
        at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26), border = brewer.pal(9,'Set1')[5:4],
        ylab='Total Variation Distance from Truth', xlab='Cells', main='Simulation Accuracy',
        axes = FALSE, ylim=c(0,max(cbind(d100,d200,d500,d1K,d2K,d5K,d10K,d20K,d50K)[1:2,])), pch=16)
axis(2)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5), labels=c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4))
abline(v=c(3,6,9,12,15,18,21,24), lty=2, col='grey50')
legend('bottomleft', fill=brewer.pal(9,'Set1')[4:5], border=NA, legend=c('Unique','EM'))



plot(c(100,5e4), range(c(0,medians[1:2,])), col='white', log = 'x',
     ylab='Median Total Variation Distance from Truth', xlab='Cells', main='Simulation Accuracy')
lines(c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4), medians[2,], col = brewer.pal(9,'Set1')[4], type='b', pch=16)
lines(c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4), medians[1,], col = brewer.pal(9,'Set1')[5], type='b', pch=16)
legend('bottomleft',col=brewer.pal(9,'Set1')[4:5], legend=c('Unique','EM'), pch=16,bty='n')

dev.off()



# Kullback-Leibler Divergence
##############################

layout(matrix(1:2,nrow=1))

boxplot(cbind(t(d100)[,3:4],t(d200)[,3:4],t(d500)[,3:4],t(d1K)[,3:4],t(d2K)[,3:4],t(d5K)[,3:4],t(d10K)[,3:4],t(d20K)[,3:4],t(d50K)[,3:4]),
        at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26), border = brewer.pal(9,'Set1')[5:4],
        ylab='Kullback-Leibler Divergence from Truth', xlab='Cells', main='Simulation Accuracy',
        axes = FALSE, ylim=c(0,max(cbind(d100,d200,d500,d1K,d2K,d5K,d10K,d20K,d50K)[3:4,])), pch=16)
axis(2)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5), labels=c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4))
abline(v=c(3,6,9,12,15,18,21,24), lty=2, col='grey50')
legend('bottomleft', fill=brewer.pal(9,'Set1')[4:5], border=NA, legend=c('Unique','EM'))



plot(c(100,5e4), range(c(0,medians[3:4,])), col='white', log = 'x',
     ylab='Median Kullback-Leibler Divergence from Truth', xlab='Cells', main='Simulation Accuracy')
lines(c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4), medians[4,], col = brewer.pal(9,'Set1')[4], type='b', pch=16)
lines(c(100,200,500,1e3,2e3,5e3,1e4,2e4,5e4), medians[3,], col = brewer.pal(9,'Set1')[5], type='b', pch=16)
legend('bottomleft',col=brewer.pal(9,'Set1')[4:5], legend=c('Unique','EM'), pch=16,bty='n')





























# old version



plot(log(c(100,20000))+c(-1,1), range(cbind(d100,d500,d1K,d5K,d10K,d20K)[1:2,]), col='white',
     ylab='Total Variation Distance from Truth', xlab='Sample Size', main='Simulation Results',
     axes = FALSE)
#abline(h=seq(.5,.8,by=.05), lty=2, col=rgb(0,0,0,.3))
axis(2)
axis(1, at=log(c(100,500,1e3,5e3,1e4,2e4)), labels = c(100,500,1e3,5e3,1e4,2e4))
boxplot(t(d100[1:2,]), add=TRUE, at=log(100)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d500[1:2,]), add=TRUE, at=log(500)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d1K[1:2,]), add=TRUE, at=log(1000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d5K[1:2,]), add=TRUE, at=log(5000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d10K[1:2,]), add=TRUE, at=log(10000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d20K[1:2,]), add=TRUE, at=log(20000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
legend('bottomleft',fill=brewer.pal(9,'Set1')[4:5], legend=c('Unique','EM'), bty='n')




plot(log(c(100,20000))+c(-1,1), range(cbind(d100,d500,d1K,d5K,d10K,d20K)[3:4,]), col='white',
     ylab='Kullback-Leibler Divergence from Truth', xlab='Sample Size', main='Simulation Results',
     axes = FALSE)
#abline(h=seq(.5,.8,by=.05), lty=2, col=rgb(0,0,0,.3))
axis(2)
axis(1, at=log(c(100,500,1e3,5e3,1e4,2e4)), labels = c(100,500,1e3,5e3,1e4,2e4))
boxplot(t(d100[3:4,]), add=TRUE, at=log(100)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d500[3:4,]), add=TRUE, at=log(500)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d1K[3:4,]), add=TRUE, at=log(1000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d5K[3:4,]), add=TRUE, at=log(5000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d10K[3:4,]), add=TRUE, at=log(10000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])
boxplot(t(d20K[3:4,]), add=TRUE, at=log(20000)+c(-.15,.15), boxwex=.3, axes = FALSE, col = brewer.pal(9,'Set1')[4:5])


legend('bottomleft',fill=brewer.pal(9,'Set1')[4:5], legend=c('Unique','EM'))


