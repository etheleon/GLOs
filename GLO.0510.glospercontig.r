#!/usr/bin/env Rscript

library(igraph)
library(arcdiagram)
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
theme_set(theme_minimal())
args = commandArgs(T)

load(sprintf("data/%s.rda",args[1]))
readglo 	= tbl_df(read.table(sprintf("~/sequencing_output/out/seq.0241/%s-family-GLOS",args[1]),comment.char="",fill=T,h=T))
contigmapping 	= tbl_df(setNames(read.table(sprintf("~/GLO/data/454readstatus/%s",args[1]),comment.char=""), c('readID', "status","startcontig","start5","startstrand","endcontig","end3","endstrand",'Overlap')))
contiglength 	= tbl_df(setNames(read.table(sprintf("data/%s-contiglength",args[1])), c("Overlap","contig","length")))
glotable = 	tbl_df(setNames(read.table(sprintf("~/sequencing_output/out/seq.0250/%s",args[1])), c("KO","GLO","reads","genuslca","family")))

#Merges GLO with Read
contigmapping = tbl_df(merge(contigmapping,readglo,by='readID'))
contigmapping$glolength = sapply(contigmapping$GLO, function(x) { length(unlist(strsplit(as.character(x),'_')))})

contigmapping2 	= 	contigmapping %.% 
    			group_by(startcontig) %.% 
    			summarise(uniqglos = length(unique(GLO)), mean_glo_length = mean(glolength))

#Plotting the number of unique GLOs against the length of the GLOs###################################
pdf(sprintf("out/glo.0510/%s_uniqueGLOs.per.contig.pdf", args[1]))
qplot(uniqglos, meanglolength , data = contigmapping2)+
	scale_x_log10()+
    	xlab("# unique GLOs")+
    	ylab("Average length of GLOs")
dev.off()
####################################################################################################

colnames(contigmapping)[colnames(contigmapping) == 'startcontig']<-'contig'
contigmapping=merge(contigmapping, contiglength, by=c("Overlap","contig"))

#The distances of the reads
m=as.matrix(dist(dm,method='binary'))
m2 <- melt(m)[melt(upper.tri(m))$value,]

trimmed_contigmapping = cbind(
select(contigmapping, c(Overlap:readID, GLO:length)), 
do.call(rbind,apply(contigmapping, 1, function(x) { 
startend=sort(c(as.integer(x[names(x) == 'start5']),as.integer(x[names(x) == 'end3'])))
    data.frame(
    start=startend[[1]], 
    end=startend[[2]]
    )
})))

trimmed_contigmapping = trimmed_contigmapping %.% 
	mutate(overlapcontig = paste(Overlap, contig,sep=""))

#The number of reads per contig
readspercontig = 	trimmed_contigmapping %.% 
	group_by(overlapcontig) %.% 
	summarise(no.of.reads = n()) %.% 
	arrange(no.of.reads)
pdf(sprintf("out/glo.0510/.pdf", args[1]))
qplot(no.of.reads, data = la4)
dev.off()


lapply(unique(trimmed_contigmapping$Overlap), function(overlap) { 
    pdf(sprintf("out/glo.0510/%s.pdf",overlap),w=10,h=5)
    lapply(unique(subset(trimmed_contigmapping, Overlap == overlap)$contig), function(contigg){ 
    singlecontig = filter(trimmed_contigmapping, Overlap == as.character(overlap), contig == as.character(contigg))


if(unique(singlecontig$length) > 200) {
rangelist = setNames(
		    data.frame(t(sapply(seq(1,as.integer(unique(singlecontig$length)), by=100), function(x) data.frame(start = x, end = x+99)))), 
		    c("start","end")
		    )


gloloc = do.call(rbind,(apply(singlecontig, 1, function(x) { 
	GLO = x[names(x) == 'GLO']
	readID = x[names(x) == 'readID']
	middle = as.integer(x[names(x) == 'start']) + ((as.integer(x[names(x) == 'end']) - as.integer(x[names(x) == 'start']))/2)
	df = 
	rangelist[apply(rangelist, 1, function(x) {
	    if(middle >= x[names(x) == 'start'] & middle <= x[names(x) == 'end'] ){ T }else{ F }
	    }),]
	    if(nrow(df) >0) { 
	df$GLO = GLO
	df$readID = readID
	df$loc = middle 
	df
	    }
	})))

glotable$num = 1:nrow(glotable)

vertex.details	= merge(gloloc, select(glotable, c(GLO,num)), by='GLO')
vertex.details$node = paste(vertex.details$start, vertex.details$end, sep="-")

#Get the relavant distances

    vertices = vertex.details %.% 
    group_by(node) %.% 
    summarise(readnum = n())

vertices = vertices[match(apply(rangelist, 1, function(x) paste(x[1],x[2],sep="-")), vertices$node),]
vertices = vertices[complete.cases(vertices),]


if(length(unique(vertex.details$node)) > 2) { 

#For all possible pairs of the nodes, calculate the distance 
edgelist = as.matrix(subset(
do.call(rbind,
apply(t(combn(unique(vertex.details$node), 2)), 1, function(x) {

#GLOs in that node
	pair1 = subset(vertex.details, node == x[1])$num
	pair2 = subset(vertex.details, node == x[2])$num

	data.frame(
	node1=x[1],
	node2=x[2],
	shared = nrow(subset(unique(melt(m[pair1,pair2])), value < 0.2))
	)
})), shared > 0))

if(nrow(edgelist) > 0) { 
g = graph.data.frame(edgelist, vertices = vertices,directed =F)

arcplot(get.edgelist(g),vertices = vertices$node, lwd.arcs = E(g)$shared, cex.nodes = sqrt(V(g)$readnum), lty.arcs = )

    }}}
    })
   dev.off() 
    })
