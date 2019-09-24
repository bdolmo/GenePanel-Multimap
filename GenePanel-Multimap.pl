#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd;
use Sort::Key::Natural qw(natsort);

my $bed;
my $readL;
my $offset;
my $mismatch;
my $genome;
my $indexed;
my $nCPU;
my $outdir;
my $devNull    = ">/dev/null 2>&1";

 # Edit line
 my $gem_mapper = `which gem-mapper`;
 chomp $gem_mapper;

if (!$gem_mapper) {
	print " ERROR: gem-mapper was not found on PATH\n";
}

my $bedtools = `which bedtools`;
chomp $bedtools;

if (!$bedtools) {
	print " ERROR: bedtools was not found on PATH\n";
}

my $Rscript = `which Rscript`;
chomp $Rscript;
if (!$Rscript) {
	print " ERROR: Rscript was not found on PATH\n";
}

 invocation () if (@ARGV <1 or !GetOptions(
  'bed:s' 	    =>\$bed,
  'l=i'		    =>\$readL, 
  'm=i'		    =>\$mismatch, 	#% mismatch allowed for mapping (default 1%)
  'g=s'		    =>\$genome, 
  'ncpu=i'      =>\$nCPU,
  'index:s' 	=>\$indexed,
  'o:s'		    =>\$outdir,
  )
 );
 if (!$outdir ) {
	print "ERROR: no output directory (-outdir) introduced\n";
	exit;
 }

 if (!-e $outdir ) {
	mkdir $outdir;
 }

 sub invocation {
 	print("
 Usage:  $0 -bed <region_of_interest> -l <read_length> -offset <bases> -m <\%mismatch> -g <genome> -index <gem_index> -temp <yes|no>\n
 Params:

  -bed    FILE   BED file of regions 
  -l      INT    Bait length (default 80)
  -m      INT    Maximum percentage of mismatches for mapping (default 1%)
  -g      FILE   Genome file in FASTA
  -ncpu   INT    Number of CPUs used for gem-mapper (default 1)
  -index  FILE   PATH to GEM indexed file
  -o      FILE   Output directory\n\n");
 exit;
 }

 $readL	     = 101  if !$readL;
 $mismatch   = 1    if !$mismatch;
 $mismatch   = 1    if $mismatch > 100;
 $nCPU       = 1    if !$nCPU;
 $offset     = 0    if !$offset; 

 print  "INFO: " . localtime() . "\n";
 print  "INFO: ROI file $bed\n";
 print  "INFO: Read length of $readL bp\n";
 print  "INFO: $mismatch % mismatches allowed\n";

 my $name = basename ($bed);
 $name =~s/\.bed//;

 my $bedmodified = "$name.modified.bed"; 
 my $bedAndFasta = "$name.fasta"; 
 my $fastq       = "$name.fastq";
 
 our %exonPos = ();
 our %summary = ();
 our %offsetInfo = ();
 our %geneHash = ();

 readBed2Struct($bed, $genome);

 # Adding padding bases at offset
 my $cmd = "cat $bed | awk '{ print \$1\"\t\"\$2-$offset\"\t\"\$3+$offset\"\t\"\$4}' > $outdir/$bedmodified";
 system ("$cmd");

 # Getting fasta sequences for every targeted region
 $cmd = "$bedtools getfasta -fi $genome -tab -name -bed $outdir/$bedmodified -fo $bedAndFasta";
 system ($cmd);

 # Writing FASTQ and extracting GC content
 open (FASTQ, ">", "$outdir/$fastq") || die "Unable to open $outdir/$fastq\n";
 open (FASTA, "<", $bedAndFasta) || die "Unable to open $bedAndFasta\n";
 while (my $line=<FASTA>) {

	my ($exon, $chr, $startOff, $endOff, $seq) = split (/[\t:-]+/, $line);

	for (my $i = 0; $i < length($seq)-$readL; $i++) {

		my $subseq = substr($seq, $i, $readL);
		my $qual = "I" x $readL;
		print FASTQ "\@$exon,$i,$chr,$startOff,$endOff\n$subseq\n+\n$qual\n";

		$offsetInfo{$exon}{$i+1}{GC} = getGC($subseq);
		$offsetInfo{$exon}{TOTAL}++;
	}
 }
 close FASTA;
 close FASTQ;

 ##########################3##################outdir########
 # Mapping baits with GEM			     			 #					
 # NOTE: important to use a multiple-hit-aware mapper# 
 # to locate all hits (including non-optimal)        #
 #####################################################

 print  "INFO: Mapping sequences with gem-mapper\n";
 
 $cmd = "$gem_mapper -I $indexed -i $outdir/$fastq -o $outdir/$name -m $mismatch -C -q offset-33 -T $nCPU $devNull";
 system($cmd) if !-e "$outdir/$name.map";

 analyzeMultimap();

 plotLowMapExons("$outdir/$name.multimap_results.txt");

 unlink($bedmodified, $bedAndFasta, $fastq);

 print "INFO: Done\n";

####### END #######

sub  plotLowMapExons {
	my $infile = shift;
	my $outfile = $infile;
	$outfile =~s/.multimap_results.txt/.lowMappability.png/;
	open (R, ">", "$outdir/plotLowMapExons.R");
	print R "library(ggplot2)\n";
	print R "mydata<-read.table(file = \"$infile\", sep = \"\t\", header=TRUE)\n";
	print R "attach(mydata)\n";
	print R "lowmap <- mydata[ which(mydata\$Mappability<100), ]\n";
	print R "myplot<-ggplot(lowmap, aes(x=reorder(lowmap\$Exon, -lowmap\$Mappability), y=lowmap\$Mappability)) + geom_bar(stat=\"identity\", fill = \"darkgreen\" ) + xlab(\"Exons\") + ylab(\"Mappability(%)\")+ theme(axis.text.x = element_text(angle = 70, hjust = 1)) + coord_flip() + theme_bw()\n";
	print R "png(\"$outfile\", res = 200, height=1400, width= 1400)\n";
	print R "myplot\n";
	print R "dev.off()\n";
	`Rscript $outdir/plotLowMapExons.R`;
	close R;
}
#######################################
#  Fill a hash with exon coordinates  #
#######################################

sub readBed2Struct {

	my $bed = shift;
	my $reference = shift;
	my $gcBed = $bed;
	$gcBed =~s/.bed/.gc.bed/;

	my $cmd = "$bedtools nuc -fi $reference -bed $bed > $gcBed";
	system $cmd;

	open (BED, "<", $gcBed) || die " ERROR: Unable to open $gcBed\n";
	while (my $line =<BED>) {
		chomp $line;
		next if $line =~/6_pct_gc/;
		my ($chr, $start, $end, $info, $at, $gc) = split /\t/, $line;
		
		$gc = sprintf "%.2f", 100*$gc;
		#NM_000256_5_6;MYBPC3
		my @tmp = split (/[;_]/, $info);
		my $gene = $tmp[-1];

		my $exonNum = $tmp[3];

		$geneHash{$gene}{$exonNum}{CHR} = $chr;
		$geneHash{$gene}{$exonNum}{START} = $start;
		$geneHash{$gene}{$exonNum}{END} = $end;
		$geneHash{$gene}{$exonNum}{GC}  = $gc;
		$exonPos{$info}{CHR}   = $chr;
		$exonPos{$info}{START} = $start;
		$exonPos{$info}{END}   = $end;
		$exonPos{$info}{GC}  = $gc;
	}
	close BED;
}

#######################
#  Analyze multimaps  #
#######################

sub analyzeMultimap {

   my $tag = 0;
   my $flag = 0;
   my $is_same_gene;

   my %tag = ();
   $tag{M} = 0;
   $tag{U} = 0;

   open (MAP, "<", "$outdir/$name.map") || die "Unable to open $outdir/$name.map\n";
   while (my $line =<MAP>) {
	chomp $line;
	my @tmp = split /\t/, $line;
	my @tmpHits = split (/,/, $tmp[4]);
   	my ($exon, $pos, $chrHead, $startHead, $endHead) = split (/,/, $tmp[0]);
   	$summary{$exon}{$pos}{HITS} = scalar @tmpHits;
   }
   close MAP;
 
	open LOG, ">", "$outdir/$name.multimap_results.txt" || die "ERROR: Unable to open $outdir/$name.multimap_results.txt\n";
	print LOG "chr\tstart\tend\tExon\tMappability\tGC\n";

	foreach my $exon (natsort keys %summary) {

		my $multiMaps = 0;
		my $uniqueMaps = 0;
		my $totalPos = 0;
		foreach my $position ( natsort keys %{$summary{$exon}} ) {
			if ($summary{$exon}{$position}{HITS} > 1) {
				$multiMaps++;
			}
			else {
				$uniqueMaps++;
			}
			$totalPos++;
		}
		my $mappability =  sprintf "%.2f", 100* $uniqueMaps/$totalPos;
		print LOG "$exonPos{$exon}{CHR}\t$exonPos{$exon}{START}\t$exonPos{$exon}{END}\t$exon\t$mappability\t$exonPos{$exon}{GC}\n";
   }
}
 ################################
 #  Plot multimapping per exon  #
 ################################

 sub plotMultimapping {

	my $infile = shift;
	my $exon   = shift;
	my $st     = shift;
	my $ends   = shift;
	my $str    = shift;
	my $genome_start = shift;
	my $genome_end = shift;

	my @tmp = split /[;_]/, $exon;
	my $exonName = $tmp[3];
	my $gene     = $tmp[4];

	open (R, ">", "$outdir/toplot.r");
	print R "library(ggplot2)\n";
	print R "mydata<-read.table(file =\"$infile\", sep = \"\t\", header = TRUE)\n";
	print R "attach(mydata)\n";
	print R "maxHits<-max(Hits)\n";
	print R "png(\"$outdir/$exon\_st\_ends.png\", res = 250, height= 1200, width=2400)\n";
	print R "myplot<-ggplot(mydata,aes(x=Relative_position, y=Hits)) + geom_line(colour=\"red\") + ylim(0, maxHits+1) + theme_minimal() + theme(panel.border = element_rect(colour = \"black\", size = 0.8, fill=NA))+ geom_vline(xintercept =0, colour=\"blue\")+ geom_vline(xintercept =$ends, colour=\"green\")+ labs(title = \"$gene Exon $exonName\") + geom_rect(data=mydata, mapping=aes(xmin=0, xmax=$ends, ymin=0.2, ymax=0.8), fill=\"dark blue\", alpha=0.2) + annotate(geom=\"text\", x=$ends/2, y=0.50, label=\"Exon $exonName\", color=\"white\") + annotate(geom=\"text\", x=0, y=maxHits+0.5, label=\"$genome_start\", color=\"black\", angle = 35) + annotate(geom=\"text\", x=$ends, y=maxHits+0.5, label=\"$genome_end\", color=\"black\", angle = 35)\n";
	print R "myplot + xlab(\"Relative position\")\n";
	print R "dev.off()\n";
	close R;

	`$Rscript $outdir/toplot.r`;
	unlink ("$outdir/toplot.r");
 }

######################################
#  Obtain GC content from a sequence #
######################################

 sub getGC {
	my $seq = shift;
	my $count = 0;
	my @tmp = split //, $seq;
	foreach my $ntd (@tmp) {
		$count++ if $ntd =~/[G|g|C|c]/;
	}
	my $perc = sprintf "%.2f", 100*$count/length($seq);
	return $perc;
 }