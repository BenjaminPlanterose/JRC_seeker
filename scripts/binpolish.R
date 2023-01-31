#Rscript scripts/binpolish.R -s {input.segments} -b {params.blacklist} -m {params.mappability} 
# -l {input.labels} -y {params.segment_min_sz} -i {params.binpolishflank} -o {output.bp} -r {output.reg}

# Read and separate arguments
options = commandArgs(trailingOnly = TRUE)
index_names = startsWith(options, '-')
arguments = options[!index_names]
names(arguments) = options[index_names]
segments = arguments['-s']
blacklist = arguments['-b']
mappability = arguments['-m']
labels = arguments['-l']
segment_min_sz = as.integer(arguments['-y'])
epsilon = as.integer(arguments['-i'])
output.bp = arguments['-o']
output.reg = arguments['-r']


############################ 1. Reading files ############################
# segments = '/media/ultron/2tb_disk1/ben/JRC_project/0_testingJRCseeker/chromhmm/output_files/test_run_4_segments.bed'
# blacklist = '/home/ultron/opt/JRC_seeker/assets/blacklist_regions/hg38_blacklist_regions.bed'
# mappability = '/home/ultron/opt/JRC_seeker/assets/mappability_files/hg38_k100.bismap.bed'
# labels = '/media/ultron/2tb_disk1/ben/JRC_project/0_testingJRCseeker/binpolish/assets/state_labels.txt'
# segment_min_sz = 200
# epsilon = 1000L
# Out
# output.bp = arguments['-o'] # polished segmentation
# output.reg = arguments['-r'] # im_regions.txt

############################ 0. Read files ############################

library(data.table)
library(GenomicRanges)
segmentations = fread(segments)
colnames(segmentations) = c('chr', 'start', 'end', 'state')
#segmentations$cpg_count = fread(cpg)$V5
#seg_counts_um = fread(um_processed)
#segmentations$M_count = seg_counts_um$V5
#segmentations$U_count = seg_counts_um$V6
labl_assign = fread(labels)
labels = labl_assign$label
names(labels) = labl_assign$state
segmentations$state = factor(segmentations$state, levels = c('E1', 'E2', 'E3', 'E4'))
levels(segmentations$state) = unname(labels[levels(segmentations$state)])
bl = fread(blacklist)
map = fread(mappability, skip = 1)

############################ 1. Build GRanges ############################

bl = reduce(GRanges(seqnames = bl$V1, ranges = IRanges(start = bl$V2, end = bl$V3)))
map = reduce(GRanges(seqnames = map$V1, ranges = IRanges(start = map$V2, end = map$V3)))

IM = GRanges(seqnames = segmentations$chr[state = segmentations$state == 'IM'], 
             ranges = IRanges(start = segmentations$start[state = segmentations$state == 'IM'], 
                              end = segmentations$end[state = segmentations$state == 'IM']))
U = GRanges(seqnames = segmentations$chr[state = segmentations$state == 'U'], 
             ranges = IRanges(start = segmentations$start[state = segmentations$state == 'U'], 
                              end = segmentations$end[state = segmentations$state == 'U']))
M = GRanges(seqnames = segmentations$chr[state = segmentations$state == 'M'], 
            ranges = IRanges(start = segmentations$start[state = segmentations$state == 'M'], 
                             end = segmentations$end[state = segmentations$state == 'M']))
ND = GRanges(seqnames = segmentations$chr[state = segmentations$state == 'ND'], 
            ranges = IRanges(start = segmentations$start[state = segmentations$state == 'ND'], 
                             end = segmentations$end[state = segmentations$state == 'ND']))
complete = reduce(Reduce(c, list(M, U, IM, ND)))

############################ 2. Extend IM regions by epsilon ############################

start(IM) = start(IM) - epsilon
end(IM) = end(IM) + epsilon
IM = reduce(IM)
U = setdiff(U, IM)
M = setdiff(M, IM)

############################ 3. Remove lop map/blacklist ############################

IM = setdiff(intersect(IM, map), bl)
IM = IM[width(IM) > segment_min_sz]
U = setdiff(intersect(U, map), bl)
M = setdiff(intersect(M, map), bl)
ND = setdiff(complete, Reduce(c, list(M, U, IM)))

############################ 4. Recombine and export ############################

IM$label = 'IM'
M$label = 'M'
U$label = 'U'
ND$label = 'ND'
polished_seg = sort(Reduce(c, list(M, IM, U, ND)))

#setwd('/media/ultron/2tb_disk1/ben/JRC_project/0_testingJRCseeker/')
polished_seg = as.data.frame(polished_seg)[,-(4:5)]
polished_seg$start = as.integer(polished_seg$start)
polished_seg$end = as.integer(polished_seg$end)

fwrite(polished_seg, file = output.bp, sep = '\t', col.names = F)

#setwd('/media/ultron/2tb_disk1/ben/JRC_project/0_testingJRCseeker/')
IM_regions = paste(seqnames(IM), paste(start(IM), end(IM), sep = '-'), sep = ':')
fwrite(data.frame(regions = IM_regions), file = output.reg, sep = '\t', col.names = F)



