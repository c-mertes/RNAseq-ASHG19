rv <- fread('/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/allelic_counts.tsv.gz')

g <- ggplot(rv, aes(refCount, altCount)) + geom_point(alpha = .7, color = 'gray70') + 
       scale_y_log10() + scale_x_log10() + theme_bw() + labs(x = 'REF count', y = 'ALT count')
g
# Add 10 counts filter
rv[, ymax := c(10:0, rep(0, .N - 11))]
rv[, x := 0:(.N - 1)]
g1 <- g + geom_ribbon(aes(ymin = 0, ymax = ymax, x = x),
                      fill = 'steelblue1', alpha = .7)
g1

# Remove Reference higher
rv[, xall:=sort(refCount)]
rv[, ymax2 := pmax(ymax, 0)]
g2 <- g1 + geom_ribbon(aes(x=xall, ymin=ymax, ymax=xall), 
                fill="steelblue", alpha = .7)
g2

# Remove non significant ones
g3 <- g2 + geom_ribbon(aes(x=xall, ymin=xall, ymax=xall+10 + xall*0.2), 
                       fill="steelblue4", alpha = .7)
g3

# Add rare
rv[, rare := F]
rv[variantID %in% c('rs369499936', 'rs3738476', 'rs759267778'), rare := T]
g4 <- g3 + geom_point(aes(col = rare), size = 2) + 
  scale_color_manual(values = c('gray70', 'firebrick')) #+ theme(legend.position = "none")
g4


library(cowplot)
plot_grid(g1, g2, g3, g4, nrow = 2)
