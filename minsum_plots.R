library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(wesanderson)
library(magrittr)
library(scales)

setwd('~/hapcut/Results')

grand <- c(wes_palette("GrandBudapest")[3],wes_palette("GrandBudapest")[2])
darjeeling <- wes_palette('Darjeeling')

#read in minsum files
cg100 <- fread('100_hapcut.parsed.txt') 
cg101 <- fread('101_hapcut.parsed.txt')
cg102 <- fread('102_hapcut.parsed.txt')
cg107 <- fread('107_hapcut.parsed.txt')
cg27 <- fread('27_hapcut.parsed.txt')
cg28 <- fread('28_hapcut.parsed.txt')
cg141 <- fread('141.hapcut.parsed.txt') 
cg193 <- fread('193.hapcut.parsed.txt')
cg33 <- fread('33.hapcut.parsed.txt')
cg108 <- fread('108.hapcut.parsed.txt')
cg198 <- fread('198.hapcut.parsed.txt') 
cg174 <- fread('174.hapcut.parsed.txt')
cg41 <- fread('41.hapcut.minsum')

#relatedness
cg196 <- fread('196.hapcut.minsum')
cg87 <- fread('87.hapcut.minsum')
cg119 <- fread('119.hapcut.minsum')
cg120 <- fread('120.hapcut.minsum')
cg145 <- fread('145.hapcut.minsum')
cg60 <- fread('60.hapcut.minsum')

#species-wide
sw103 <- fread('103_hapcut.minsum') 
sw5a <- fread('5a_hapcut.minsum')
sw83 <- fread('83_hapcut.minsum')
sw85 <- fread('85_hapcut.minsum')
sw86 <- fread('86_hapcut.minsum')
sw88 <- fread('88_hapcut.minsum')
sw91 <- fread('91_hapcut.minsum')
sw93 <- fread('93_hapcut.minsum')
sw94 <- fread('94_hapcut.minsum')
sw95 <- fread('95_hapcut.minsum')
sw97 <- fread('97_hapcut.minsum')
swP8 <- fread('P8_hapcut.minsum')
swaxe <- fread('axe_hapcut.minsum')


#filter for locus
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
cg41.FT <- cg41 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)

cg119.FT <- cg119 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg120.FT <- cg120 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg145.FT <- cg145 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg196.FT <- cg196 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg60.FT <- cg60 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
cg87.FT <- cg87 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)

sw103.FT <- sw103 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw5a.FT <- sw5a %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw83.FT <- sw83 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw85.FT <- sw85 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw86.FT <- sw86 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw88.FT <- sw88 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw91.FT <- sw91 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw93.FT <- sw93 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw94.FT <- sw94 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw95.FT <- sw95 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
sw97.FT <- sw97 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
swaxe.FT <- swaxe %>% filter(chr=='scaffold_2', last>7539046, first<7541032)
swP8.FT <- swP8 %>% filter(chr=='scaffold_2', last>7539046, first<7541032)


#read in data that spans the locus
cg100.scent <- fread('../Data//100.scent')
cg101.scent <- fread('../Data/101.scent')
cg102.scent <- fread('../Data/102.scent')
cg107.scent <- fread('../Data/107.scent')
cg108.scent <- fread('../Data/108.scent')
cg141.scent <- fread('../Data/141.scent')
cg174.scent <- fread('../Data/174.scent')
cg193.scent <- fread('../Data/193.scent')
cg198.scent <- fread('../Data/198.scent')
cg27.scent <- fread('../Data/27.scent')
cg28.scent <- fread('../Data/28.scent')
cg33.scent <- fread('../Data/33.scent')
cg41.scent <- fread('../Data/41.scent')
cg87.scent <- fread('../Data/87.scent')
cg196.scent <- fread('../Data/196.scent')

#plot minsum per block
cg100.plot <- cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '100') + geom_point(data = cg100.FT, aes(x=phased, y=minsum), color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg100.pdf",cg100.plot)
cg101.plot <- cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '101') + geom_point(data = cg101.FT, inherit.aes=TRUE, color=grand[2])  + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg101.pdf",cg101.plot)
cg102.plot <- cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '102') + geom_point(data = cg102.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg102.pdf",cg102.plot)
cg107.plot <- cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '107') + geom_point(data = cg107.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg107.pdf",cg107.plot)
cg196.plot <- cg196 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '196') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg196.pdf",cg196.plot)
cg27.plot <- cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '27') + geom_point(data = cg27.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg27.pdf",cg27.plot)
cg28.plot <- cg28 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '28') + geom_point(data = cg28.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg28.pdf",cg28.plot)
cg87.plot <- cg87 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '87') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg87.pdf",cg87.plot)
cg198.plot <- cg198 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '198') + geom_point(data = cg198.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg198.pdf",cg198.plot)
cg141.plot <- cg141 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '141') + geom_point(data = cg141.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg141.pdf",cg141.plot)
cg193.plot <- cg193 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '193') + geom_point(data = cg193.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg193.pdf",cg193.plot)
cg33.plot <- cg33 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '33') + geom_point(data = cg33.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg33.pdf",cg33.plot)
cg108.plot <- cg108 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '108') + geom_point(data = cg108.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg108.pdf",cg108.plot)
cg174.plot <- cg174 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '174') + geom_point(data = cg174.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg174.pdf",cg174.plot)
cg41.plot <- cg41 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '41') + geom_point(data = cg41.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg41.pdf",cg41.plot)
cg87.plot <- cg87 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '87') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg87.pdf",cg87.plot)

#panel plot of the above
cg100.plot <- cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '100') + geom_point(data = cg100.FT, aes(x=phased, y=minsum), color=grand[2]) + xlab('') + 
  ylab('')

cg101.plot <- cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '101') + geom_point(data = cg101.FT, inherit.aes=TRUE, color=grand[2])  + xlab('') + 
  ylab('')

cg102.plot <- cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '102') + geom_point(data = cg102.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('minimum pairwise distance to C. rubella')

cg107.plot <- cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '107') + geom_point(data = cg107.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('')

cg196.plot <- cg196 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '196') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('')

cg27.plot <- cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '27') + geom_point(data = cg27.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('')

cg141.plot <- cg141 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '141') + geom_point(data = cg141.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('')

cg193.plot <- cg193 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '193') + geom_point(data = cg193.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('')

cg33.plot <- cg33 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '33') + geom_point(data = cg33.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('')

cg108.plot <- cg108 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '108') + geom_point(data = cg108.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('minimum pairwise distance to C. rubella')

cg174.plot <- cg174 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '174') + geom_point(data = cg174.FT, inherit.aes=TRUE, color=grand[2]) + xlab('') + 
  ylab('')

panel.1 <- plot_grid(cg100.plot, cg101.plot, cg102.plot, cg107.plot, cg27.plot, labels=c('A','','','','',''),ncol=1,nrow=5,align='v')
panel.2 <- plot_grid(cg33.plot, cg141.plot, cg108.plot, cg174.plot, cg193.plot, labels=c('B','','','','',''),ncol=1,nrow=5,align='v')

panel.plot <- plot_grid(panel.1,panel.2,ncol=2)
save_plot("panel.pdf",panel.plot,nrow=5,ncol=2,base_aspect_ratio=1.3)

#plot minsum per block and within locus
cg100.both.plot <- cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '100') + geom_point(data = cg100.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(47,23),color=darjeeling[3])
save_plot("cg100both.pdf",cg100.both.plot)
cg101.both.plot <- cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '101') + geom_point(data = cg101.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(39,12),color=darjeeling[3])
save_plot("cg101both.pdf",cg101.both.plot)
cg102.both.plot <- cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '102') + geom_point(data = cg102.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(37,11),color=darjeeling[3])
save_plot("cg102both.pdf",cg102.both.plot)
cg107.both.plot <- cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '107') + geom_point(data = cg107.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(31,8),color=darjeeling[3])
save_plot("cg107both.pdf",cg107.both.plot)
cg108.both.plot <- cg108 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '108') + geom_point(data = cg108.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(39,0),color=darjeeling[3])
save_plot("cg108both.pdf",cg108.both.plot)
cg141.both.plot <- cg141 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '141') + geom_point(data = cg141.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(34,0),color=darjeeling[3])
save_plot("cg141both.pdf",cg141.both.plot)
cg174.both.plot <- cg174 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '174') + geom_point(data = cg174.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(34,0),color=darjeeling[3])
save_plot("cg174both.pdf",cg174.both.plot)
cg193.both.plot <- cg193 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '193') + geom_point(data = cg193.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(33,0),color=darjeeling[3])
save_plot("cg193both.pdf",cg193.both.plot)
cg198.both.plot <- cg198 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '198') + geom_point(data = cg198.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(47,22),color=darjeeling[3])
save_plot("cg198both.pdf",cg198.both.plot)
cg27.both.plot <- cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '27') + geom_point(data = cg27.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(35,13),color=darjeeling[3])
save_plot("cg27both.pdf",cg27.both.plot)
cg33.both.plot <- cg33 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '33') + geom_point(data = cg33.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(35,0),color=darjeeling[3])
save_plot("cg33both.pdf",cg33.both.plot)
cg41.both.plot <- cg41 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '41') + geom_point(data = cg41.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella') +geom_point(aes(31,10),color=darjeeling[3])
save_plot("cg41both.pdf",cg41.both.plot)

#make panel plot
#plot minsum per block and within locus
cg100.both.plot <- cg100 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '100') + geom_point(data = cg100.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(47,23),color=darjeeling[3])

cg101.both.plot <- cg101 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '101') + geom_point(data = cg101.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(39,12),color=darjeeling[3])

cg102.both.plot <- cg102 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '102') + geom_point(data = cg102.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('minimum pairwise distance to C. rubella') +geom_point(aes(37,11),color=darjeeling[3])

cg107.both.plot <- cg107 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '107') + geom_point(data = cg107.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(31,8),color=darjeeling[3])

cg108.both.plot <- cg108 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '108') + geom_point(data = cg108.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('minimum pairwise distance to C. rubella') +geom_point(aes(39,0),color=darjeeling[3])

cg141.both.plot <- cg141 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '141') + geom_point(data = cg141.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(34,0),color=darjeeling[3])

cg174.both.plot <- cg174 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '174') + geom_point(data = cg174.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(34,0),color=darjeeling[3])

cg193.both.plot <- cg193 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '193') + geom_point(data = cg193.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('') +geom_point(aes(33,0),color=darjeeling[3])

cg198.both.plot <- cg198 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '198') + geom_point(data = cg198.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(47,22),color=darjeeling[3])

cg27.both.plot <- cg27 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '27') + geom_point(data = cg27.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('phased SNPS/block') + 
  ylab('') +geom_point(aes(35,13),color=darjeeling[3])

cg33.both.plot <- cg33 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = darjeeling[2]) +
  ggtitle(label = '33') + geom_point(data = cg33.FT, aes(x=phased, y=minsum), color=darjeeling[1]) + xlab('') + 
  ylab('') +geom_point(aes(35,0),color=darjeeling[3])

panel.both.1 <- plot_grid(cg100.both.plot, cg101.both.plot, cg102.both.plot, cg107.both.plot, cg27.both.plot, labels=c('A','','','','',''),ncol=1,nrow=5,align='v')
panel.both.2 <- plot_grid(cg33.both.plot, cg141.both.plot, cg108.both.plot, cg174.both.plot, cg193.both.plot, labels=c('B','','','','',''),ncol=1,nrow=5,align='v')

panel.both.plot <- plot_grid(panel.both.1,panel.both.2,ncol=2)
save_plot("panel.both.pdf",panel.both.plot,nrow=5,ncol=2,base_aspect_ratio=1.3)

#read in simulation files
conserv5 <- fread('gr.sim.conservative.5.zerosum')
conserv20 <- fread('gr.sim.conservative.20.zerosum')
conserv100 <- fread('gr.sim.conservative.100.zerosum')
btwn5 <- fread('gr.sim.conservative.5.minsumdist')
btwn20 <- fread('gr.sim.conservative.20.minsumdist')
btwn100 <- fread('gr.sim.conservative.100.minsumdist')
within5 <- fread('gr.sim.conservative.5.withindist')
within20 <- fread('gr.sim.conservative.20.withindist')
within100 <- fread('gr.sim.conservative.100.withindist')

#minsum plots
btwn5.hist <- ggplot(btwn5, aes(minsum, fill=grand[1])) +
  geom_histogram(bins=80)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('pairwise distance to C. rubella')
save_plot('btwn5.hist.pdf',btwn5.hist)

btwn20.hist <- ggplot(btwn20, aes(minsum, fill=grand[1])) +
  geom_histogram(bins=125)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('pairwise distance to C. rubella')
save_plot('btwn20.hist.pdf',btwn20.hist)

btwn100.hist <- ggplot(btwn100, aes(minsum, fill=grand[1])) +
  geom_histogram(bins=99)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('pairwise distance to C. rubella')
save_plot('btwn100.hist.pdf',btwn100.hist)

#zero divergence pdf plots
conserv5.hist <- ggplot(conserv5, aes(zeros, fill=grand[1]))+
  geom_histogram(bins=30)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('haplotypes in simulation with no divergence\n from simulated reference')
save_plot('conserv5.hist.pdf',conserv5.hist)

conserv20.hist <- ggplot(conserv20, aes(zeros, fill=grand[1]))+
  geom_histogram(bins=30)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('haplotypes in simulation with no divergence\n from simulated reference')
save_plot('conserv20.hist.pdf',conserv20.hist)

conserv100.hist <- ggplot(conserv100, aes(zeros, fill=grand[1]))+
  geom_histogram()+
  theme(legend.position='none')+
  xlab('haplotypes in simulation with no divergence\n from simulated reference')
save_plot('conserv100.hist.pdf',conserv100.hist)

#withinsum plots
within5.hist <- ggplot(within5, aes(withinsum, fill=grand[1])) +
  geom_histogram(bins=77)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('Within population minsum')+
  xlab('pairwise distance to C. grandiflora')
save_plot('within5.hist.pdf',within5.hist)

within20.hist <- ggplot(within20, aes(withinsum, fill=grand[1])) +
  geom_histogram(bins=125)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('Within population minsum')+
  xlab('pairwise distance to C. grandiflora')
save_plot('within20.hist.pdf',within20.hist)

within100.hist <- ggplot(within100, aes(withinsum, fill=grand[1])) +
  geom_histogram(bins=100)+
  theme(legend.position='none')+
  scale_y_continuous(labels=comma)+
  xlab('Within population minsum')+
  xlab('pairwise distance to C. grandiflora')
save_plot('within100.hist.pdf',within100.hist)

#final plots
minsum.plot <- plot_grid(cg100.plot, cg193.plot, labels = c("A", "B"), ncol = 2, align = "h")
save_plot('final.minsum.pdf',minsum.plot, base_aspect_ratio=1.3,
          ncol = 2)
total.plot <- plot_grid(minsum.plot, btwn100.hist, labels = c('','C'), ncol =1 , nrow=2)
save_plot('total.plot.pdf',total.plot,nrow=2,ncol=2,base_aspect_ratio=1.3)

## relatedness plots
cg119.plot <- cg119 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '119') + geom_point(data = cg119.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg119.pdf",cg119.plot)
cg120.plot <- cg120 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '120') + geom_point(data = cg120.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg120.pdf",cg120.plot)
cg145.plot <- cg145 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '145') + geom_point(data = cg145.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg145.pdf",cg145.plot)
cg196.plot <- cg196 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '196') + geom_point(data = cg196.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg196.pdf",cg196.plot)
cg60.plot <- cg60 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '60') + geom_point(data = cg60.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg60.pdf",cg60.plot)
cg87.plot <- cg87 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '87') + geom_point(data = cg87.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("cg87.pdf",cg87.plot)

## species-wide plots
sw103.plot <- sw103 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '103') + geom_point(data = sw103.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw103.pdf",sw103.plot)
sw5a.plot <- sw5a %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '5a') + geom_point(data = sw5a.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw5a.pdf",sw5a.plot)
sw83.plot <- sw83 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '83') + geom_point(data = sw83.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw83.pdf",sw83.plot)
sw85.plot <- sw85 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '85') + geom_point(data = sw85.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw85.pdf",sw85.plot)
sw86.plot <- sw86 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '86') + geom_point(data = sw86.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw86.pdf",sw86.plot)
sw88.plot <- sw88 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '88') + geom_point(data = sw88.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw88.pdf",sw88.plot)
sw91.plot <- sw91 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '91') + geom_point(data = sw91.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw91.pdf",sw91.plot)
sw93.plot <- sw93 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '93') + geom_point(data = sw93.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw93.pdf",sw93.plot)
sw94.plot <- sw94 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '94') + geom_point(data = sw94.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw94.pdf",sw94.plot)
sw95.plot <- sw95 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '95') + geom_point(data = sw95.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw95.pdf",sw95.plot)
sw97.plot <- sw97 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = '97') + geom_point(data = sw97.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("sw97.pdf",sw97.plot)
swaxe.plot <- swaxe %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = 'axe') + geom_point(data = swaxe.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("swaxe.pdf",swaxe.plot)
swP8.plot <- swP8 %>% ggplot(aes(x=phased, y=minsum)) + geom_point() + geom_smooth(method='lm', se = TRUE, color = grand[1]) +
  ggtitle(label = 'P8') + geom_point(data = swP8.FT, inherit.aes=TRUE, color=grand[2]) + xlab('phased SNPS/block') + 
  ylab('minimum pairwise distance\n to C. rubella')
save_plot("swP8.pdf",swP8.plot)
