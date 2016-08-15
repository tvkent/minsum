library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(wesanderson)

setwd('~/hapcut/Results')

grand <- c(wes_palette("GrandBudapest")[3],wes_palette("GrandBudapest")[2])

cg100 <- fread('100_hapcut.parsed.txt') 
cg101 <- fread('101_hapcut.parsed.txt')
cg102 <- fread('102_hapcut.parsed.txt')
cg107 <- fread('107_hapcut.parsed.txt')
cg196 <- fread('196_hapcut.parsed.txt')
cg27 <- fread('27_hapcut.parsed.txt')
cg28 <- fread('28_hapcut.parsed.txt')
cg87 <- fread('87_hapcut.parsed.txt')
cg141 <- fread('141.hapcut.parsed.txt') %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))
cg193 <- fread('193.hapcut.parsed.txt') %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))
cg33 <- fread('33.hapcut.parsed.txt') %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))
cg108 <- fread('108.hapcut.parsed.txt')  %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))
cg198 <- fread('198.hapcut.parsed.txt') %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))
cg174 <- fread('174.hapcut.parsed.txt') %>% filter(chr=='scaffold_2') %>% mutate(mid=((last+first)/2))

cg100.FT <- cg100 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg101.FT <- cg101 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg102.FT <- cg102 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg107.FT <- cg107 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg196.FT <- cg196 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg27.FT <- cg27 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg28.FT <- cg28 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg87.FT <- cg87 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)

cg141.FT <- cg141 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg193.FT <- cg193 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg33.FT <- cg33 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg174.FT <- cg174 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg198.FT <- cg198 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg108.FT <- cg108 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)

#get average minsum in FS locus
hapcut174 <- fread('~/hapcut/Data/174.formatted.hapcut',sep = '\t')
hapcut33 <- fread('~/hapcut/Data/33.formatted.hapcut',sep = '\t')
hapcut141 <- fread('~/hapcut/Data/141.formatted.hapcut',sep = '\t')
hapcut193 <- fread('~/hapcut/Data/193.formatted.hapcut',sep = '\t')
hapcut108 <- fread('~/hapcut/Data/108.formatted.hapcut',sep = '\t')

hapcut174.filter <- hapcut174 %>% filter(V4=='scafold_2')
FS.min.174 <- hapcut174 %>% filter(V4=='scafold_2',V5>=7539046 & V5<=7541032)# %>% summarize(sum1=sum(V2),sum2=sum(V3)) %>% min()
hapcut33.filter <- hapcut33 %>% filter(V4=='scafold_2')
FS.min.33 <- hapcut33 %>% filter(V4=='scafold_2',V5>=7539046 & V5<=7541032)# %>% summarize(sum1=sum(V2),sum2=sum(V3)) %>% min()
hapcut141.filter <- hapcut141 %>% filter(V4=='scafold_2')
FS.min.141 <- hapcut141 %>% filter(V4=='scafold_2',V5>=7539046 & V5<=7541032)# %>% summarize(sum1=sum(V2),sum2=sum(V3)) %>% min()
hapcut193.filter <- hapcut193 %>% filter(V4=='scafold_2')
FS.min.193 <- hapcut193 %>% filter(V4=='scafold_2',V5>=7539046 & V5<=7541032)# %>% summarize(sum1=sum(V2),sum2=sum(V3)) %>% min()
hapcut141.filter <- hapcut141 %>% filter(V4=='scafold_2')
FS.min.108 <- hapcut108 %>% filter(V4=='scafold_2',V5>=7539046 & V5<=7541032)# %>% summarize(sum1=sum(V2),sum2=sum(V3)) %>% min()

#plot minsum along ch2
pdf('minsum.scaffold2.pdf',width=8, height=5)
plot(cg174$minsum~cg174$mid, pch=19)
lines(smooth.spline(x = cg174$mid, y=cg174$minsum,spar = .3),col='red')
points(y=FS.min.174,x=7540039,pch=19,col='blue')
plot(cg141$minsum~cg141$mid, pch=19)
lines(smooth.spline(x = cg141$mid, y=cg141$minsum,spar = .3),col='red')
points(y=FS.min.141,x=7540039,pch=19,col='blue')
plot(cg193$minsum~cg193$mid, pch=19)
lines(smooth.spline(x = cg193$mid, y=cg193$minsum,spar = .3),col='red')
points(y=FS.min.193,x=7540039,pch=19,col='blue')
plot(cg33$minsum~cg33$mid, pch=19)
lines(smooth.spline(x = cg33$mid, y=cg33$minsum,spar = .3),col='red')
points(y=FS.min.33,x=7540039,pch=19,col='blue')
plot(cg108$minsum~cg108$mid, pch=19)
lines(smooth.spline(x = cg108$mid, y=cg108$minsum,spar = .3),col='red')
points(y=FS.min.108,x=7540039,pch=19,col='blue')
dev.off()

#what percent of blocks of FT block size are same or smaller?
ratio.141 <- ((cg141 %>% filter(phased==50) %>% nrow())/(cg141 %>% filter(phased==50, minsum<=0)
                                                         %>% nrow()))/100
ratio.193 <- ((cg193 %>% filter(phased==41) %>% nrow())/(cg193 %>% filter(phased==41, minsum<=0)
                                                         %>% nrow()))/100
ratio.33 <- ((cg33 %>% filter(phased==101) %>% nrow())/(cg33 %>% filter(phased==101, minsum<=3)
                                                         %>% nrow()))/100
ratio.174 <- ((cg174 %>% filter(phased==139) %>% nrow())/(cg174 %>% filter(phased==139, minsum<=15)
                                                         %>% nrow()))/100
ratio.108 <- ((cg108 %>% filter(phased==58) %>% nrow())/(cg108 %>% filter(phased==58, minsum<=8)
                                                        %>% nrow()))/100

#all <- rbind(cg100, cg101, cg102, cg107, cg196, cg27, cg28, cg87)

#plot span by minsum
# pdf('haps.pdf',width=8, height=5)
# cg100 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '100') + geom_point(data = cg100.FT, inherit.aes=TRUE, color=grand[2])
# cg101 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '101') + geom_point(data = cg101.FT, inherit.aes=TRUE, color=grand[2])
# cg102 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '102') + geom_point(data = cg102.FT, inherit.aes=TRUE, color=grand[2])
# cg107 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '107') + geom_point(data = cg107.FT, inherit.aes=TRUE, color=grand[2])
# cg196 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '196') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2])
# cg27 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '27') + geom_point(data = cg27.FT, inherit.aes=TRUE, color=grand[2])
# cg28 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '28') + geom_point(data = cg28.FT, inherit.aes=TRUE, color=grand[2])
# cg87 %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '87') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2])
# all %>% ggplot(aes(x=span, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = 'All') + geom_point(data = all.FT, inherit.aes=TRUE, color=grand[2])
# dev.off()
# 
# #plot SNPs by minsum
# pdf('haps_snps.pdf',width=8, height=5)
# cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '100') + geom_point(data = cg100.FT, inherit.aes=TRUE, color=grand[2])
# cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '101') + geom_point(data = cg101.FT, inherit.aes=TRUE, color=grand[2])
# cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '102') + geom_point(data = cg102.FT, inherit.aes=TRUE, color=grand[2])
# cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '107') + geom_point(data = cg107.FT, inherit.aes=TRUE, color=grand[2])
# cg196 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '196') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2])
# cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '27') + geom_point(data = cg27.FT, inherit.aes=TRUE, color=grand[2])
# cg28 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '28') + geom_point(data = cg28.FT, inherit.aes=TRUE, color=grand[2])
# cg87 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = '87') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2])
# all %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
#   ggtitle(label = 'All') + geom_point(data = all.FT, inherit.aes=TRUE, color=grand[2])
# dev.off()

pdf('FT_minsum.pdf',width=8, height=5)
cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '100-NI') + geom_point(data = cg100.FT, inherit.aes=TRUE, color=grand[2])
cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '101-NI') + geom_point(data = cg101.FT, inherit.aes=TRUE, color=grand[2])
cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '102-NI') + geom_point(data = cg102.FT, inherit.aes=TRUE, color=grand[2])
cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '107-NI') + geom_point(data = cg107.FT, inherit.aes=TRUE, color=grand[2])
cg196 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '196-NI') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2])
cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '27-NI') + geom_point(data = cg27.FT, inherit.aes=TRUE, color=grand[2])
cg28 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '28-NI') + geom_point(data = cg28.FT, inherit.aes=TRUE, color=grand[2])
cg87 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '87-NI') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2])
cg198 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '198-NI-green') + geom_point(data = cg198.FT, inherit.aes=TRUE, color=grand[2])
cg141 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '141') + geom_point(data = cg141.FT, inherit.aes=TRUE, color=grand[2])
cg193 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '193') + geom_point(data = cg193.FT, inherit.aes=TRUE, color=grand[2])
cg33 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '33') + geom_point(data = cg33.FT, inherit.aes=TRUE, color=grand[2])
cg108 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '108') + geom_point(data = cg108.FT, inherit.aes=TRUE, color=grand[2])
cg174 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '174') + geom_point(data = cg174.FT, inherit.aes=TRUE, color=grand[2])
dev.off()
