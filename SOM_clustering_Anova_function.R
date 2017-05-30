som.cluster <- function(xdim = 4, ydim = 4, topo = c("rectangular", "hexagonal"), transformation = 'asin', method = c('aov', 'lm', 'fisher'), plot = TRUE, size.aj = TRUE, manual = FALSE){
  
  path.read <- getwd()
  if(dir.exists(paste0(path.read, '/output'))){
    path.store <- paste0(path.read, '/output/')
  } else {
    dir.create(paste0(path.read, '/output'))
    path.store <- paste0(path.read, '/output/')
  }
  # List of packages for session
  .packages = c("ggplot2", "kohonen","dplyr", "tidyr", "NMF", "colorspace")
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
  
  # Load packages into session 
  lapply(.packages, require, character.only=TRUE)
  
  #Load files automatically in the workspace
  temp <- list.files(pattern = '*.csv')
  myfiles <- lapply(temp, read.csv)
  
  #Add groups using file names, perform transformation, downsample if needed
  
  if(size.aj == TRUE){ #This if loop is if samples have to be checked for size adjustment
    nE <- unlist(lapply(myfiles, FUN = nrow))
    down.sample <- which(nE > median(nE)*1.5)
    for(i in down.sample){
      myfiles[[i]] <- myfiles[[i]][sample(x = nrow(myfiles[[i]]),size = ceiling(median(nE))),]
    }
    #Check if transformation is desired and if yes perform it
    if(transformation == 'asin'){
      cof <- readline('What coefficient you want to use for arcsin transformation [standard cof = 5]? ')
      cof <- strtoi(cof)
      if(cof == ''){
        cof <- 5
      }
      for(i in 1:length(myfiles)){
        myfiles[[i]][myfiles[[i]] == 0] <- runif(length(which(myfiles[[i]] == 0)), min = -1, max = 0)
        myfiles[[i]] <- asinh(myfiles[[i]]/cof)
        myfiles[[i]]$Group1 <- strsplit(temp[[i]], "_")[[1]][[1]]
        myfiles[[i]]$Group2 <- strsplit(temp[[i]], "_")[[1]][[2]]
        myfiles[[i]]$Group3 <- paste0('Sample_', strsplit(strsplit(temp[[i]], "_")[[1]][[3]], "\\.")[[1]][[1]])
        myfiles[[i]]$Group4 <- i
      }
    } else { #does not perform transformation
      for(i in 1:length(myfiles)){
        myfiles[[i]][myfiles[[i]] == 0] <- runif(length(which(myfiles[[i]] == 0)), min = -1, max = 0)
        myfiles[[i]]$Group1 <- strsplit(temp[[i]], "_")[[1]][[1]]
        myfiles[[i]]$Group2 <- strsplit(temp[[i]], "_")[[1]][[2]]
        myfiles[[i]]$Group3 <- paste0('Sample_', strsplit(strsplit(temp[[i]], "_")[[1]][[3]], "\\.")[[1]][[1]])
        myfiles[[i]]$Group4 <- i
      }
    }
  } else{  #Adjustment of size not chekced and performed
    if(transformation == 'asin'){ # transformation step if selected
      cof <- readline('What coefficient you want to use for arcsin transformation [standard cof = 5]? ')
      cof <- strtoi(cof)
      if(cof == ''){
        cof <- 5
      }
      for(i in 1:length(myfiles)){
        myfiles[[i]][myfiles[[i]] == 0] <- runif(length(which(myfiles[[i]] == 0)), min = -1, max = 0)
        myfiles[[i]] <- asinh(myfiles[[i]]/cof)
        myfiles[[i]]$Group1 <- strsplit(temp[[i]], "_")[[1]][[1]]
        myfiles[[i]]$Group2 <- strsplit(temp[[i]], "_")[[1]][[2]]
        myfiles[[i]]$Group3 <- paste0('Sample_', strsplit(strsplit(temp[[i]], "_")[[1]][[3]], "\\.")[[1]][[1]])
        myfiles[[i]]$Group4 <- i
      }
    } else { # no transformation performed
      for(i in 1:length(myfiles)){
        myfiles[[i]][myfiles[[i]] == 0] <- runif(length(which(myfiles[[i]] == 0)), min = -1, max = 0)
        myfiles[[i]]$Group1 <- strsplit(temp[[i]], "_")[[1]][[1]]
        myfiles[[i]]$Group2 <- strsplit(temp[[i]], "_")[[1]][[2]]
        myfiles[[i]]$Group3 <- paste0('Sample_', strsplit(strsplit(temp[[i]], "_")[[1]][[3]], "\\.")[[1]][[1]])
        myfiles[[i]]$Group4 <- i
      }
    }
    }
  
  # Create a data frame for further analysis
  df <- do.call(rbind, myfiles)
  
  #Select markers for SOM analysis
  columns <- list()
  print('Select markers you want to use for SOM analysis.')
  
  for(i in 1:ncol(df)){
    columns[[i]] <- readline(paste0('Load ', colnames(df)[i], '?', ' Y/n '))
  }
  columns <- do.call(rbind, columns)
  
  #Data for SOM
  df_analysis <- df[which(columns == "y" | columns == "Y" |columns == "")]
  
  #Keep grouping variables separately, not really used later 
  grouping <- df[grep('Group', colnames(df))]
  
  #Perform SOM analysis
  set.seed(27) #for reporoducibility
  som.clas <- kohonen::som(as.matrix(df_analysis), grid = somgrid(xdim,ydim,topo))
  
  #Add clusters to original data 
  df <- cbind(df, Cluster=som.clas$unit.classif)
  
  #Grouping by groups and Clusters, and counting how many events are found in each subgroup
  df_s <- dplyr::group_by(df, Group1, Group2, Group3, Cluster) %>% summarise(n_events=n())
  
  #Table with each sample in separate column
  df_s1 <- tidyr::spread(df_s, Group3, n_events)
  
  group1 <- readline('Type which grouping (1 or 2) variable you want to use [most often 2]: ')

  #####################
  #
  #Stopped here with commenting 
  df_s1[is.na(df_s1)] <- 0
  if(nrow(df_s1) != ydim*as.numeric(count(unique(df_s1[paste0("Group",group1)])))){
    #do magic
    #group_by(df_s1, Group2) %>% summarise(n())
    for(i in 1:as.numeric(count(unique(df_s1[paste0("Group",group1)])))){
    df_temp <- df_s1[eval(parse(text=paste0("which(df_s1$Group",group1," == unique(df_s1$Group",group1,")[i])"))),]
      if(nrow(df_temp) < xdim*ydim){
      temp1 <- as.data.frame(matrix(nrow = xdim*ydim - nrow(df_temp), ncol = ncol(df_temp)))
      colnames(temp1) <- colnames(df_temp)
        for(a in 1:length(which(1:(xdim*ydim) %in% df_temp$Cluster == FALSE))){
        temp1[a,1:2] <- df_temp[1, 1:2]
        temp1[a,3:ncol(temp1)] <- 0
        temp1[a, "Cluster"] <- which(1:(xdim*ydim) %in% df_temp$Cluster == FALSE)[a]
        df_tem1 <- as.data.frame(df_s1)
        df_s1 <- as.data.frame(rbind(df_tem1, temp1))
        }
      } # end of if statement first
    }
  }# end of if statement second
  
  df_s2 <- tidyr::gather(data = df_s1,"Group3", "n_events", grep("Samp", colnames(df_s1)))
  df_s2 <- df_s2[c(colnames(df_s2)[grep("Group", colnames(df_s2))], colnames(df_s2)[-grep("Group", colnames(df_s2))])]
  if (length(method) > 1){
   #write code that will work if multiple methods are selected e.g. anova and fisher
    for(i in 1:length(method)){
      if(method[[i]] == 'lm' | method[[i]] == 'aov'){
      x <- unique(df_s2$Cluster)
      models <- sapply(x, function(x) {
        bartlett.test(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster == x)
      }, simplify=FALSE)
      
      hom <- cbind(Cluster = x, Bartlet.p = sapply(models, function(x) {x$p.value}))
      }
      if(method[[i]] == 'lm'){
        x <- unique(df_s2$Cluster)
        models <- sapply(x, function(x) {
          lm(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster==x)
        }, simplify=FALSE)
      }
      if(method[[i]] == 'aov'){
        x <- unique(df_s2$Cluster)
        models <- sapply(x, function(x) {
          aov(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster == x)
        }, simplify=FALSE)
      }
      if(method[[i]] == 'fisher'){
        #x <- unique(df_s2$Cluster)
        fisher <- sapply(split(df_s1, f = df_s1$Cluster), function(x) {
          if(nrow(x) >=2){
            m <- as.matrix(x[grep('Sample', colnames(x))])
            fisher.test(x = m, hybrid = T)$p.value
          }
        }, simplify=FALSE)
      }
      if(method[[i]] == 'aov' | method[[i]] == 'lm'){
        ANOVA.tables <- sapply(models, anova, simplify=FALSE)
        df_f1 <- do.call(rbind, ANOVA.tables)
        
        df_f1[grep('Group', row.names(df_f1)),ncol(df_f1)+1] <- hom[,2]
        colnames(df_f1)[ncol(df_f1)] <- 'Bartlett.p'
        df_f1[grep('Group', row.names(df_f1)),ncol(df_f1)+1] <- x
        colnames(df_f1)[ncol(df_f1)] <- 'Cluster'
      }
      if(method[[i]] == 'fisher'){
        fisher <- do.call(rbind, fisher)
        df_f2 <- data.frame(Cluster = unique(df_s2$Cluster), Fisher.p = fisher)
      }
    }
  }else{
  if(method == 'lm' | method == 'aov'){
  x <- unique(df_s2$Cluster)
  models <- sapply(x, function(x) {
    bartlett.test(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster == x)
  }, simplify=FALSE)

  hom <- cbind(Cluster = x, Bartlet.p = sapply(models, function(x) {x$p.value}))
  }

  if(method == 'lm'){
    x <- unique(df_s2$Cluster)
    models <- sapply(x, function(x) {
      lm(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster==x)
    }, simplify=FALSE)
  }
  if(method == 'aov'){
    x <- unique(df_s2$Cluster)
    models <- sapply(x, function(x) {
      aov(as.formula(paste0('n_events ~ Group', group1)), data=df_s2, subset = Cluster == x)
    }, simplify=FALSE)
  }
  if(method == 'fisher'){
    #x <- unique(df_s2$Cluster)
    fisher <- sapply(split(df_s1, f = df_s1$Cluster), function(x) {
      if(nrow(x) >=2){
        m <- as.matrix(x[grep('Sample', colnames(x))])
        fisher.test(x = m, hybrid = T)$p.value
      }
    }, simplify=FALSE)
  }
  
  if(method == 'aov' | method == 'lm'){
  ANOVA.tables <- sapply(models, anova, simplify=FALSE)
  df_f <- do.call(rbind, ANOVA.tables)
  
  df_f[grep('Group', row.names(df_f)),ncol(df_f)+1] <- hom[,2]
  colnames(df_f)[ncol(df_f)] <- 'Bartlett.p'
  df_f[grep('Group', row.names(df_f)),ncol(df_f)+1] <- x
  colnames(df_f)[ncol(df_f)] <- 'Cluster'
  }
  if(method == 'fisher'){
  fisher <- do.call(rbind, fisher)
  df_f <- data.frame(Cluster = unique(df_s2$Cluster), Fisher.p = fisher)
  }
  }
   
  #Create plots for all clusters 
  #if(plot == TRUE && max(x) <= 25){
  p <-ggplot(df_s2, aes_string(paste0('Group', group1), 'n_events', fill=paste0('Group', group1)))+
      facet_grid(. ~ Cluster)+
      geom_boxplot()+
      geom_point()+
      theme_classic()
  
  pdf(paste0(path.store, 'Plot1.pdf'), w = length(unique(df_s2$Cluster))*1, h=5)
  print(p)
  dev.off()
  #} 
  
  # df:       transformed original data with SOM clusters appended at the end
  # hom:      homogeneity, Bartlett's test, p has to be > 0.05 for ANOVA to be valid
  # df_f:     ANOVA or Fisher test summary table for all clusters
  # df_s:     Summary table with number of events (data points) in all groups and clusters
  # df_s1:    Summary table with number of events (data points) in all groups and clusters, seen for each sample
  # som.clas: SOM output... description later
  # models:   List of all linear models or aov models 
  # plot:     ggplot form ready to be plotted
  
  if (exists('df_f1') & exists('df_f2')){
    total <- list(df_fA=df_f1, df_fF=df_f2, df_s = df_s2, df_s1 = df_s1, df=df, som.clas=som.clas, models=models, plot=p)
    if(manual == T){
      namefile <- list()
      namefile[[1]] <- readline('Type in name for ANOVA test summary table: ')
      namefile[[2]] <- readline('Type in name for Fisher test summary table: ')
      namefile[[3]] <- readline('Type in name for summary table, number of events in each cluster split to groups: ')
      namefile[[4]] <- readline('Type in name for summary table, number of events in each cluster split to groups and samples: ')
      namefile[[5]] <- readline('Type in name for transformed original data with SOM clusters appended at the end: ')
      for(i in c(1,2,3,4,5)){
        if(i == 1){
          write.csv(total[[i]], paste0(path.store, namefile[[i]], "_Table_", i,".csv"), row.names = T)
        }
        write.csv(total[[i]], paste0(path.store, namefile[[i]], "_Table_", i,".csv"), row.names = F)
      }
    } else {
      for(i in c(1,2,3,4,5)){
        if(i == 1){
          write.csv(total[[i]], paste0(path.store,"Table_", i,".csv"), row.names = T)
        }
        write.csv(total[[i]], paste0(path.store,"Table_", i,".csv"), row.names = F)
      }
    }  
    } else {
  total <- list(df_f=df_f, df_s = df_s2, df_s1 = df_s1, df=df, som.clas=som.clas, models=models, plot=p)
  if(manual == T){
    namefile <- list()
    namefile[[1]] <- readline('Type in name for ANOVA or Fisher test summary table: ')
    namefile[[2]] <- readline('Type in name for summary table, number of events in each cluster split to groups: ')
    namefile[[3]] <- readline('Type in name for summary table, number of events in each cluster split to groups and samples: ')
    namefile[[4]] <- readline('Type in name for transformed original data with SOM clusters appended at the end: ')
    for(i in c(1,2,3,4)){
      if(i == 1){
        write.csv(total[[i]], paste0(path.store, namefile[[i]], "_Table_", i,".csv"), row.names = T)
      }
      write.csv(total[[i]], paste0(path.store, namefile[[i]], "_Table_", i,".csv"), row.names = F)
    }
  } else {
    for(i in c(1,2,3,4)){
      if(i == 1){
        write.csv(total[[i]], paste0(path.store, "Table_", i,".csv"), row.names = T)
      }
      write.csv(total[[i]], paste0(path.store, "Table_", i,".csv"), row.names = F)
    }
  }
  }
  
  heatm <- readline('Do you want to generate heatmap of clusters and groups [Y/n]? ')
  if (heatm == 'Y'){
    columns_p <- list()
    print('Select markers you want to use for SOM analysis.')
    
    for(i in 1:ncol(df)){
      columns_p[[i]] <- readline(paste0('Load ', colnames(df)[i], '?', ' Y/n '))
    }
    columns_p <- do.call(rbind, columns_p)
  }
    
  df_sum <- aggregate(df[which(columns_p == "y" | columns_p == "Y" |columns_p == "")], list(df$Group2, df$Cluster), median)
  write.csv(df_sum, paste0(path.store, "Average_markers_in_clusters_groups.csv"), row.names = F)
  
  iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100, 
                       plot=FALSE, random=FALSE) {
    # Presently doesn't allow hmax > hmin (H is circular)
    # n: number of colours
    # hmin: lower bound of hue (0-360)
    # hmax: upper bound of hue (0-360)
    # cmin: lower bound of chroma (0-180)
    # cmax: upper bound of chroma (0-180)
    # lmin: lower bound of luminance (0-100)
    # lmax: upper bound of luminance (0-100)
    # plot: plot a colour swatch?
    # random: should clustering be random? (if FALSE, seed will be set to 1,
    #         and the RNG state will be restored on exit.) 
    require(colorspace)
    stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, 
              hmax <= 360, cmax <= 180, lmax <= 100, 
              hmin <= hmax, cmin <= cmax, lmin <= lmax,
              n > 0)
    if(!random) {
      if (exists(".Random.seed", .GlobalEnv)) {
        old_seed <- .GlobalEnv$.Random.seed
        on.exit(.GlobalEnv$.Random.seed <- old_seed)
      } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
      }
      set.seed(1)
    }
    lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), 
                                     seq(-100, 100, 5), 
                                     seq(-110, 100, 5))))
    if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
             hmax != 360 || cmax != 180 || lmax != 100))) {
      hcl <- as(lab, 'polarLUV')
      hcl_coords <- coords(hcl)
      hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                         hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin & 
                         hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
      #hcl <- hcl[-which(is.na(coords(hcl)[, 2]))]
      lab <- as(hcl, 'LAB')    
    }
    lab <- lab[which(!is.na(hex(lab))), ]
    clus <- kmeans(coords(lab), n, iter.max=50)
    if (isTRUE(plot)) {
      swatch(hex(LAB(clus$centers)))
    }
    hex(LAB(clus$centers))
  }
  
  cluster_colors <- iwanthue(xdim*ydim)
  
  #Var1 = c("navy", "darkgreen")
  names(cluster_colors) = as.factor(x)
  Var2 = iwanthue(length(unique(eval(parse(text=paste0("df_s1$Group",group1))))))
  names(Var2) = as.factor(unique(eval(parse(text=paste0("df_s1$Group",group1)))))
  ann_colors = list(Group1 = Var2, Group2 = cluster_colors)
  #aheatmap(x, annCol = annotation, annColors = ann_colors)
  
  pdf(paste0(path.store, 'HeatPlot.pdf'))
  NMF::aheatmap(df_sum[-grep('Group', colnames(df_sum))], annRow = data.frame(Group1 = as.factor(t(df_sum$Group.1)), Group2 = as.factor(t(df_sum$Group.2))), annColors = ann_colors, fontsize=7, distfun = "pearson", cellwidth = length(columns_p)*0.2)
  dev.off()  
  
  #Writing the files in output directory
  
  # if(manual == T){
  #   namefile <- list()
  # namefile[[1]] <- readline('Type in name for ANOVA or Fisher test summary table: ')
  # namefile[[2]] <- readline('Type in name for summary table, number of events in each cluster split to groups: ')
  # namefile[[3]] <- readline('Type in name for summary table, number of events in each cluster split to groups and samples: ')
  # namefile[[4]] <- readline('Type in name for transformed original data with SOM clusters appended at the end: ')
  #   for(i in c(1,2,3,4)){
  #     if(i == 1){
  #       write.csv(total[[i]], paste0(path.store, '/', method, namefile[[i]], "_Table_", i,".csv"), row.names = T)
  #     }
  #     write.csv(total[[i]], paste0(path.store, '/', method, namefile[[i]], "_Table_", i,".csv"), row.names = F)
  #   }
  # } else {
  #   for(i in c(1,2,3,4)){
  #     if(i == 1){
  #       write.csv(total[[i]], paste0(method, "_Table_", i,".csv"), row.names = T)
  #     }
  #     write.csv(total[[i]], paste0(method, "_Table_", i,".csv"), row.names = F)
  #   }
  # }
  #Function returns all calculations in a list 'total'
  return(total)
}
