library(Rsamtools)



seq_contig <-'' 
pos <- 0

for (i in 1:length(bam_data)) {
	current_pos <- bam_data[i, 'pos']
	current_seq <- bam_data[i, 'seq']
	if (nchar(seq_contig) == 0) {
		seq_contig <- current_seq
		pos <- current_pos
		seq_length <- nchar(current_seq) 
	} else {
		pos_dif <- current_pos - pos
		if (pos_dif <= seq_length) {
			overlap_length <- pos + seq_length - current_pos
			next_seq <- substr(current_seq, overlap_length+1, nchar(current_seq))
			seq_contig <- paste0(seq_contig, next_seq)
			seq_length <- nchar(current_seq)
			pos <- current_pos
		} else {
			inter_char <- paste0('-',pos + seq_length, 'to', current_pos-1 , '-' )
			seq_contig <- paste0(seq_contig, inter_char, current_seq)
			pos <- current_pos
			seq_length <- nchar(current_seq)
		}
	}
}
