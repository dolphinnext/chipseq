$HOSTNAME = ""
params.outdir = 'results'  

// enable required indexes to build them
params.use_Bowtie2_Index = (params.run_Sequential_Mapping == "yes" || params.run_Bowtie2 == "yes") ? "yes" : ""
params.use_Bowtie_Index  = (params.run_Sequential_Mapping == "yes") ? "yes" : ""
params.use_STAR_Index    = (params.run_Sequential_Mapping == "yes") ? "yes" : ""
if (!params.reads){params.reads = ""} 
if (!params.run_Motif_Finder_on_ChIP_Peaks){params.run_Motif_Finder_on_ChIP_Peaks = ""} 
if (!params.mate){params.mate = ""} 
if (!params.genome_url){params.genome_url = ""} 
if (!params.gtf_url){params.gtf_url = ""} 
if (!params.commondb_url){params.commondb_url = ""} 

Channel
	.fromFilePairs( params.reads , size: (params.mate != "pair") ? 1 : 2 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_1_reads_g123_3;g_1_reads_g123_18}

Channel.value(params.run_Motif_Finder_on_ChIP_Peaks).set{g_119_run_Homer_g132_7}
Channel.value(params.mate).into{g_122_mate_g_68;g_122_mate_g123_3;g_122_mate_g123_11;g_122_mate_g123_16;g_122_mate_g123_18;g_122_mate_g123_19;g_122_mate_g123_20;g_122_mate_g123_21;g_122_mate_g124_26;g_122_mate_g124_30;g_122_mate_g124_32;g_122_mate_g127_10;g_122_mate_g127_13;g_122_mate_g128_9;g_122_mate_g128_23;g_122_mate_g128_25;g_122_mate_g126_82;g_122_mate_g126_95;g_122_mate_g126_123;g_122_mate_g126_126}
g_129_genome_url_g125_15 = file(params.genome_url, type: 'any') 
g_130_gtf_url_g125_15 = file(params.gtf_url, type: 'any') 
Channel.value(params.commondb_url).set{g_131_commondb_url_g125_15}

params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_1_reads_g123_18.into{g123_18_reads_g123_19}
g123_18_log_file_g123_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_1_reads_g123_18
 val mate from g_122_mate_g123_18

output:
 set val(name), file("reads/*")  into g123_18_reads_g123_19
 file "*.{fastx,trimmomatic}.log"  into g123_18_log_file_g123_11

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
system("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);



system("!{runGzip}");
my $quality="33";
my $format="";
my ($format, $len)=getFormat("!{file1}");
print "fastq format: $format\\n";
if ($format eq "solexa"){   
    $quality="64";
}
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        system("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        print "trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log";
        system("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log";
        system("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}


# automatic format detection
sub getFormat
{
   my ($filename)=@_;

   # set function variables
   open (IN, $filename);
   my $j=1;
   my $qualities="";
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $qualities.=$line;
        $len=length($line);
     }
     $j++;
   }
   close(IN);
  
   my $format = "";

   # set regular expressions
   my $sanger_regexp = qr/[!"#$%&'()*+,-.\\/0123456789:]/;
   my $solexa_regexp = qr/[\\;<=>\\?]/;
   my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\\[\\]\\^\\_\\`abcdefgh]/;
   my $all_regexp = qr/[\\@ABCDEFGHI]/;

   # set counters
   my $sanger_counter = 0;
   my $solexa_counter = 0;
   my $solill_counter = 0;

   # check qualities
   if( $qualities =~ m/$sanger_regexp/ ){
          $sanger_counter = 1;
   }
   if( $qualities =~ m/$solexa_regexp/ ){
          $solexa_counter = 1;
   }
   if( $qualities =~ m/$solill_regexp/ ){
          $solill_counter = 1;
   }

   # determine format
   if( $sanger_counter ){
        $format = "sanger";
    }elsif($solexa_counter ){
        $format = "solexa";
    }elsif($solill_counter ){
        $format = "illumina";
    }

    # return file format
    return( $format,$len );
}





'''

}
}


params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g123_18_reads_g123_19.into{g123_19_reads_g123_20}
g123_19_log_file_g123_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g123_18_reads_g123_19
 val mate from g_122_mate_g123_19

output:
 set val(name), file("reads/*")  into g123_19_reads_g123_20
 file "*.log" optional true  into g123_19_log_file_g123_21

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
system("mkdir reads");

system("!{runGzip}");
if ("!{mate}" eq "pair") {
    my $file1 = "!{file1}";
    my $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my ($format, $len)=getFormat($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $format, $len);
    my ($format, $len)=getFormat($file2);
    print "length of $file2: $len\\n";
    trimFiles($file2, $trim2, $format, $len);
} else {
    my $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my ($format, $len)=getFormat($file1);
    print "length of file1 $len\\n";
    trimFiles($file1, $trim1, $format, $len);
}

sub trimFiles
{
  my ($file, $trim, $format, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="";
    #if ($format eq "sanger")
    #{   
      $quality="-Q33";
    #}
    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      print "parameters for $file: $param \\n ";
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "$com\\n";
      if ($com eq ""){
          print "trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          system("$com");
          print "trimmer executed for $file \\n";
      }
    }

    
}

# automatic format detection
sub getFormat
{
   my ($filename)=@_;

   # set function variables
   open (IN, $filename);
   my $j=1;
   my $qualities="";
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $qualities.=$line;
        $len=length($line);
     }
     $j++;
   }
   close(IN);
  
   my $format = "";

   # set regular expressions
   my $sanger_regexp = qr/[!"#$%&'()*+,-.\\/0123456789:]/;
   my $solexa_regexp = qr/[\\;<=>\\?]/;
   my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\\[\\]\\^\\_\\`abcdefgh]/;
   my $all_regexp = qr/[\\@ABCDEFGHI]/;

   # set counters
   my $sanger_counter = 0;
   my $solexa_counter = 0;
   my $solill_counter = 0;

   # check qualities
   if( $qualities =~ m/$sanger_regexp/ ){
          $sanger_counter = 1;
   }
   if( $qualities =~ m/$solexa_regexp/ ){
          $solexa_counter = 1;
   }
   if( $qualities =~ m/$solill_regexp/ ){
          $solill_counter = 1;
   }

   # determine format
   if( $sanger_counter ){
        $format = "sanger";
    }elsif($solexa_counter ){
        $format = "solexa";
    }elsif($solill_counter ){
        $format = "illumina";
    }

    # return file format
    return( $format,$len );
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g123_19_log_file_g123_21.collect()
 val mate from g_122_mate_g123_21

output:
 file "trimmer_summary.tsv"  into g123_21_outputFileTSV_g_115

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g123_19_reads_g123_20.into{g123_20_reads_g124_32}
g123_20_log_file_g123_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g123_19_reads_g123_20
 val mate from g_122_mate_g123_20

output:
 set val(name), file("reads/*")  into g123_20_reads_g124_32
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g123_20_log_file_g123_16

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen

// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent
    
    
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
system("mkdir reads unpaired");


system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";
my $format=getFormat("!{file1}");
print "fastq format: $format\\n";
my $quality="33";
if ($format eq "sanger"){   
    $quality="33";
} elsif ($format eq "illumina" || $format eq "solexa"){
    $quality="64";
} 
print "fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        system("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        system("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        system("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
    
}


# automatic format detection
sub getFormat
{
   my ($filename)=@_;

   # set function variables
   open (IN, $filename);
   my $j=1;
   my $qualities="";
   while(my $line=<IN>)
   {
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $qualities.=$line;
     }
     $j++;
   }
   close(IN);
  
   my $format = "";

   # set regular expressions
   my $sanger_regexp = qr/[!"#$%&'()*+,-.\\/0123456789:]/;
   my $solexa_regexp = qr/[\\;<=>\\?]/;
   my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\\[\\]\\^\\_\\`abcdefgh]/;
   my $all_regexp = qr/[\\@ABCDEFGHI]/;

   # set counters
   my $sanger_counter = 0;
   my $solexa_counter = 0;
   my $solill_counter = 0;

   # check qualities
   if( $qualities =~ m/$sanger_regexp/ ){
          $sanger_counter = 1;
   }
   if( $qualities =~ m/$solexa_regexp/ ){
          $solexa_counter = 1;
   }
   if( $qualities =~ m/$solill_regexp/ ){
          $solill_counter = 1;
   }

   # determine format
   if( $sanger_counter ){
        $format = "sanger";
    }elsif($solexa_counter ){
        $format = "solexa";
    }elsif($solill_counter ){
        $format = "illumina";
    }

    # return file format
    return( $format );
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g123_20_log_file_g123_16.collect()
 val mate from g_122_mate_g123_16

output:
 file "quality_filter_summary.tsv"  into g123_16_outputFileTSV_g_115

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /adapter_removal_detailed_summary.tsv$/) "adapter_removal_detailed_summary/$filename"
}

input:
 file logfile from g123_18_log_file_g123_11.collect()
 val mate from g_122_mate_g123_11

output:
 file "adapter_removal_summary.tsv"  into g123_11_outputFileTSV_g_115
 file "adapter_removal_detailed_summary.tsv" optional true  into g123_11_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

params.gtf =  ""  //* @input
params.genome =  ""  //* @input
params.commondb =  ""  //* @input

process Check_and_Build_Module_Check_Genome_GTF {

input:
 file fasta from g_129_genome_url_g125_15
 file downGtf from g_130_gtf_url_g125_15
 val commondb_url from g_131_commondb_url_g125_15

output:
 val "${params.genome}"  into g125_15_genomePath_g125_0, g125_15_genomePath_g125_6, g125_15_genomePath_g125_8, g125_15_genomePath_g125_10, g125_15_genomePath_g125_5, g125_15_genomePath_g125_13
 val "${params.gtf}"  into g125_15_gtfPath_g125_0, g125_15_gtfPath_g125_6, g125_15_gtfPath_g125_8, g125_15_gtfPath_g125_10, g125_15_gtfPath_g125_4, g125_15_gtfPath_g125_13
 val "${params.commondb}"  into g125_15_commondb_path_g125_18

when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
genome_dir  = params.genome.substring(0, params.genome.lastIndexOf('/')) 
slashCount = commondb_url.count("/")
cutDir = slashCount - 3;

"""
downGenomePath=\$(realpath $fasta)
downGtfPath=\$(realpath $downGtf)
if [ ! -e "${params.genome}" ] ; then
    echo "${params.genome} not found"
    mkdir -p ${genome_dir}
    cp -n \$downGenomePath ${params.genome}
fi
if [ ! -e "${params.gtf}" ] ; then
    echo "${params.gtf} not found"
    mkdir -p ${gtf_dir}
    cp -n \$downGtfPath ${params.gtf}
fi
if [ ! -e "${params.commondb}" ] ; then
    echo "${params.commondb} not found"
    mkdir -p ${params.commondb}
    wget -l inf -nc -nH --cut-dirs=$cutDir -R 'index.html*' -r --no-parent --directory-prefix=${params.commondb} $commondb_url
fi

"""




}


process Check_and_Build_Module_Bowtie_Index {

input:
 val genome from g125_15_genomePath_g125_13
 val gtf from g125_15_gtfPath_g125_13

output:
 val resultDir  into g125_13_genomeIndexPath_g125_18

when:
(params.use_Bowtie_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
bowtie_build_parameters = params.Check_and_Build_Module_Bowtie_Index.bowtie_build_parameters
basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "BowtieIndex"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.2.ebwt" ] ; then
    echo "${resultDir}/${basename}.rev.2.ebwt Bowtie index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie-build ${bowtie_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir
fi
"""




}

params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_GTF2BED12 {

input:
 val gtf from g125_15_gtfPath_g125_4


when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
"""
if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${params.bed}
fi
"""




}

params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 val genome from g125_15_genomePath_g125_5


when:
params.run_checkAndBuild == "yes"

script:
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
basename_and_path  = genome.substring(0, genome.lastIndexOf('.'))

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$0,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${basename_and_path}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${basename_and_path}.chrom.sizes
fi
"""




}


process Check_and_Build_Module_Check_Build_Rsem_Index {

input:
 val genome from g125_15_genomePath_g125_10
 val gtf from g125_15_gtfPath_g125_10

output:
 val cmdAr  into g125_10_command_g125_11
 val resultDirAr  into g125_10_path_g125_11

when:
(params.use_RSEM_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
create_bowtie_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie_rsem_index
create_bowtie2_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie2_rsem_index
create_star_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_star_rsem_index
transcript_to_gene_map = params.Check_and_Build_Module_Check_Build_Rsem_Index.transcript_to_gene_map
RSEM_build_parameters = params.Check_and_Build_Module_Check_Build_Rsem_Index.RSEM_build_parameters

genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
newDirNameAr = []
cmdAr = []
resultDirAr = []
if (create_bowtie_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie') }
if (create_bowtie2_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie2') }
if (create_star_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_STAR') }

transcript_to_gene_mapText = ""
if (transcript_to_gene_map?.trim()){
    transcript_to_gene_mapText = "--transcript-to-gene-map " + transcript_to_gene_map
}

for (i = 0; i < newDirNameAr.size(); i++) {
    resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirNameAr[i]
    tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirNameAr[i]
    resultDirAr.push(resultDir)
    cmd = ""
    indexType = ""
    if (newDirNameAr[i] == 'RSEM_ref_Bowtie'){
        indexType = "--bowtie "
        checkFile = "${basenameGenome}.rev.2.ebwt" 
    } else if (newDirNameAr[i] == 'RSEM_ref_Bowtie2'){
        indexType = "--bowtie2 "
        checkFile = "${basenameGenome}.rev.2.bt2" 
    } else if (newDirNameAr[i] == 'RSEM_ref_STAR'){
        indexType = "--star "
        checkFile = "genomeParameters.txt" 
    }
    cmd = "if [ ! -e \"${resultDir}/${checkFile}\" ] ; then rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir && rsem-prepare-reference ${RSEM_build_parameters} --gtf ${gtf} ${transcript_to_gene_mapText} ${indexType} ${genome} ${basenameGenome} && cd .. && mv $tmpResultDir $resultDir; fi"
    cmdAr.push(cmd)
}


"""

"""

}


process Check_and_Build_Module_Build_Index_RSEM_run {

input:
 val resultDir from g125_10_path_g125_11.flatten()
 val command from g125_10_command_g125_11.flatten()

output:
 val resultDir  into g125_11_genomeIndexPath

script:
"""    
$command
"""
}


process Check_and_Build_Module_Hisat2_Index {

input:
 val genome from g125_15_genomePath_g125_8
 val gtf from g125_15_gtfPath_g125_8

output:
 val resultDir  into g125_8_genomeIndexPath

when:
(params.use_Hisat2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
hisat2_build_parameters = params.Check_and_Build_Module_Hisat2_Index.hisat2_build_parameters
genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
basenameGTF = gtf.substring(gtf.lastIndexOf('/')+1,gtf.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Hisat2Index"
resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

extract_splice_sites = "hisat2_extract_splice_sites.py ${gtf} > ${tmpResultDir}/${basenameGTF}.hisat2_splice_sites.txt"
extract_exons = "hisat2_extract_exons.py ${gtf}> ${tmpResultDir}/${basenameGTF}.hisat2_exons.txt"
ss = "--ss ${basenameGTF}.hisat2_splice_sites.txt"
exon = "--exon ${basenameGTF}.hisat2_exons.txt"

"""
if [ ! -e "${resultDir}/${basenameGenome}.8.ht2" ] ; then
    echo "${resultDir}/${basenameGenome}.8.ht2 Hisat2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir 
    $extract_splice_sites
    $extract_exons
    hisat2-build ${hisat2_build_parameters} $ss $exon ${genome} ${basenameGenome}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""

}

bowtie2_build_parameters = params.Check_and_Build_Module_Bowtie2_Index.bowtie2_build_parameters

process Check_and_Build_Module_Bowtie2_Index {

input:
 val genome from g125_15_genomePath_g125_6
 val gtf from g125_15_gtfPath_g125_6

output:
 val resultDir  into g125_6_genomeIndexPath_g125_18, g125_6_genomeIndexPath_g127_13

when:
(params.use_Bowtie2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:

basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Bowtie2Index"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.1.bt2" ] ; then
    echo "${resultDir}/${basename}.rev.1.bt2 Bowtie2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie2-build ${bowtie2_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""



}


process Check_and_Build_Module_STAR_Index_Check_Build {

input:
 val genome from g125_15_genomePath_g125_0
 val gtf from g125_15_gtfPath_g125_0

output:
 val resultDir  into g125_0_genomeIndexPath_g125_18

when:
(params.use_STAR_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
star_build_parameters = params.Check_and_Build_Module_STAR_Index_Check_Build.star_build_parameters
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
indexbasedir  = gtf_dir.substring(0, gtf_dir.lastIndexOf('/'))
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "STARIndex" 
resultDir = indexbasedir +"/"+ newDirName 
tmpResultDir = indexbasedir +"/_tmp_"+ newDirName
"""
if [ ! -e "${resultDir}/SA" ] ; then
    echo "STAR index not found"
    rm -rf $tmpResultDir ${resultDir} && mkdir -p $tmpResultDir && cd $tmpResultDir
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $tmpResultDir --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
    cd .. && mv $tmpResultDir $resultDir && cd ${resultDir}
    ln -s ../../main/${filename} ${filename}
    ln -s ../../main/${filename}.fai ${filename}.fai
fi
"""





}

params.gtf =  ""  //* @input
params.genome =  ""  //* @input
params.commondb =  ""  //* @input
if (!(params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes")){
g125_15_commondb_path_g125_18.into{g125_18_commondb_path_g124_32}
} else {


process Check_and_Build_Module_Check_Sequential_Mapping_Indexes {

input:
 val commondb from g125_15_commondb_path_g125_18
 val bowtieIndex from g125_13_genomeIndexPath_g125_18
 val bowtie2Index from g125_6_genomeIndexPath_g125_18
 val starIndex from g125_0_genomeIndexPath_g125_18

output:
 val commondb  into g125_18_commondb_path_g124_32

when:
params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes"

script:
"""
"""
}
}


g125_18_commondb_path_g124_32= g125_18_commondb_path_g124_32.ifEmpty(file('commondb_path', type: 'any')) 

params.run_Sequential_Mapping =   "yes"   //* @dropdown @options:"yes","no" @show_settings:"Sequential_Mapping"
params.bowtieInd_rRNA =  ""  //* @input
params.bowtieInd_ercc =  ""  //* @input
params.bowtieInd_miRNA =  ""  //* @input
params.bowtieInd_tRNA =  ""  //* @input
params.bowtieInd_piRNA =  ""  //* @input
params.bowtieInd_snRNA =  ""  //* @input
params.bowtieInd_rmsk =  ""  //* @input
params.bowtie_index =  ""  //* @input
params.bowtie2_index =  ""  //* @input
params.star_index =  ""  //* @input

//both bowtie and bowtie2 indexes located in same path
bowtieIndexes = [rRNA: params.bowtieInd_rRNA, 
                 ercc: params.bowtieInd_ercc,
                 miRNA: params.bowtieInd_miRNA,
                 tRNA: params.bowtieInd_tRNA,
                 piRNA: params.bowtieInd_piRNA,
                 snRNA: params.bowtieInd_snRNA,
                 rmsk: params.bowtieInd_rmsk]
                 
genomeIndexes = [bowtie: params.bowtie_index,
                 bowtie2: params.bowtie2_index,
                 STAR: params.star_index+"/genome"]


//_nucleicAcidType="dna" should be defined in the autofill section of pipeline header in case dna is used.
remove_duplicates = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates
remove_duplicates_based_on_UMI_after_mapping = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates_based_on_UMI_after_mapping


_select_sequence = params.Sequential_Mapping_Module_Sequential_Mapping._select_sequence
index_directory = params.Sequential_Mapping_Module_Sequential_Mapping.index_directory
name_of_the_index_file = params.Sequential_Mapping_Module_Sequential_Mapping.name_of_the_index_file
_aligner = params.Sequential_Mapping_Module_Sequential_Mapping._aligner
aligner_Parameters = params.Sequential_Mapping_Module_Sequential_Mapping.aligner_Parameters
description = params.Sequential_Mapping_Module_Sequential_Mapping.description
filter_Out = params.Sequential_Mapping_Module_Sequential_Mapping.filter_Out

desc_all=[]
description.eachWithIndex() {param,i -> 
    if (param.isEmpty()){
        desc_all[i] = name_of_the_index_file[i]
    }  else {
        desc_all[i] = param.replaceAll("[ |.|;]", "_")
    }
}
custom_index=[]
index_directory.eachWithIndex() {param,i -> 
    if (_select_sequence[i] == "genome"){
        custom_index[i] = genomeIndexes[_aligner[i]]
    }else if (_select_sequence[i] == "custom"){
        custom_index[i] = param+"/"+name_of_the_index_file[i]
    }else {
        custom_index[i] = bowtieIndexes[_select_sequence[i]]
    }
}

mapList = []
paramList = []
alignerList = []
filterList = []
indexList = []

//concat default mapping and custom mapping
mapList = (desc_all) 
paramList = (aligner_Parameters)
alignerList = (_aligner)
filterList = (filter_Out)
indexList = (custom_index)

mappingList = mapList.join(" ") // convert into space separated format in order to use in bash for loop
paramsList = paramList.join(",") // convert into comma separated format in order to use in as array in bash
alignersList = alignerList.join(",") 
filtersList = filterList.join(",") 
indexesList = indexList.join(",") 
//* @style @condition:{remove_duplicates="yes",remove_duplicates_based_on_UMI_after_mapping},{remove_duplicates="no"},{_select_sequence="custom", index_directory,name_of_the_index_file,description,_aligner,aligner_Parameters,filter_Out},{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"),_aligner,aligner_Parameters,filter_Out}  @array:{_select_sequence,_select_sequence, index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out,description} @multicolumn:{_select_sequence,_select_sequence,index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out, description},{remove_duplicates,remove_duplicates_based_on_UMI_after_mapping}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2500
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
}
//* platform
//* autofill
if (!(params.run_Sequential_Mapping == "yes")){
g123_20_reads_g124_32.into{g124_32_reads_g_68}
g124_32_bowfiles_g124_26 = Channel.empty()
g124_32_bowfiles_g_79 = Channel.empty()
g124_32_bam_file_g124_23 = Channel.empty()
g124_32_bam_file_g124_27 = Channel.empty()
g124_32_bam_index_g124_23 = Channel.empty()
g124_32_bam_index_g124_27 = Channel.empty()
g124_32_filter_g124_26 = Channel.empty()
g124_32_log_file_g124_30 = Channel.empty()
} else {


process Sequential_Mapping_Module_Sequential_Mapping {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*\/.*_sorted.bam$/) "sequential_mapping/$filename"
	else if (filename =~ /.*\/.*_sorted.bam.bai$/) "sequential_mapping/$filename"
	else if (filename =~ /.*\/.*_duplicates_stats.log$/) "sequential_mapping/$filename"
}

input:
 set val(name), file(reads) from g123_20_reads_g124_32
 val mate from g_122_mate_g124_32
 val commondb_path from g125_18_commondb_path_g124_32

output:
 set val(name), file("final_reads/*")  into g124_32_reads_g_68
 set val(name), file("bowfiles/*") optional true  into g124_32_bowfiles_g124_26, g124_32_bowfiles_g_79
 file "*/*_sorted.bam" optional true  into g124_32_bam_file_g124_23
 file "*/*_sorted.bam.bai" optional true  into g124_32_bam_index_g124_23
 val filtersList  into g124_32_filter_g124_26
 file "*/*_sorted.dedup.bam" optional true  into g124_32_bam_file_g124_27
 file "*/*_sorted.dedup.bam.bai" optional true  into g124_32_bam_index_g124_27
 file "*/*_duplicates_stats.log" optional true  into g124_32_log_file_g124_30

when:
params.run_Sequential_Mapping == "yes"

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

"""
#!/bin/bash
mkdir reads final_reads bowfiles 
if [ -n "${mappingList}" ]; then
    $runGzip
    #rename files to standart format
    if [ "${mate}" == "pair" ]; then
        mv $file1 ${name}.1.fastq 2>/dev/null
        mv $file2 ${name}.2.fastq 2>/dev/null
        mv ${name}.1.fastq ${name}.2.fastq reads/.
    else
        mv $file1 ${name}.fastq 2>/dev/null
        mv ${name}.fastq reads/.
    fi
    #sequential mapping
    k=0
    prev="reads"
    IFS=',' read -r -a paramsListAr <<< "${paramsList}" #create comma separated array 
    IFS=',' read -r -a filtersListAr <<< "${filtersList}"
    IFS=',' read -r -a indexesListAr <<< "${indexesList}"
    IFS=',' read -r -a alignersListAr <<< "${alignersList}"
    wrkDir=\$(pwd)
    for rna_set in ${mappingList}
    do
        ((k++))
        printf -v k2 "%02d" "\$k" #turn into two digit format
        mkdir -p \${rna_set}/unmapped
        cd \$rna_set
        if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
            ln -s \${wrkDir}/\${prev}/* .
            prev=\${rna_set}/unmapped
        else
            ln -s \${wrkDir}/\${prev}/* .
        fi
        genomeDir=`dirname "\${indexesListAr[\$k-1]}"`
        echo "INFO: genomeDir: \$genomeDir"
        if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -o  -e "\${indexesListAr[\$k-1]}.fa"  -o  -e "\${indexesListAr[\$k-1]}.fasta"  -o  -e "\$genomeDir/SAindex" ]; then
            if [ -e "\${indexesListAr[\$k-1]}.fa" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fa
            elif [ -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fasta
            fi
            echo "INFO: fasta: \$fasta"
            if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -a "\${alignersListAr[\$k-1]}" == "bowtie2" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.bt2 Bowtie2 index found."
            elif [ -e "\${indexesListAr[\$k-1]}.1.ebwt" -a "\${alignersListAr[\$k-1]}" == "bowtie" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.ebwt Bowtie index found."
            elif [ -e "\$genomeDir/SAindex" -a "\${alignersListAr[\$k-1]}" == "STAR" ] ; then
                echo "INFO: \$genomeDir/SAindex STAR index found."
            elif [ -e "\${indexesListAr[\$k-1]}.fa" -o  -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2-build \$fasta \${indexesListAr[\$k-1]}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    if [ -e "\${indexesListAr[\$k-1]}.gtf" ]; then
                        STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles \$fasta --sjdbGTFfile \${indexesListAr[\$k-1]}.gtf --genomeSAindexNbases 5
                    else
                        echo "WARNING: \${indexesListAr[\$k-1]}.gtf not found. STAR index is not generated."
                    fi
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie-build \$fasta \${indexesListAr[\$k-1]}
                fi
            fi
                
            if [ "${mate}" == "pair" ]; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un-conc unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq --al-conc ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam > \${k2}_${name}.bow_\${rna_set}  2>&1
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.1.fastq ${name}.2.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.1.fastq
                    mv ${name}.starUnmapped.out.mate2 unmapped/${name}.unmapped.2.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}   \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq -S  \${rna_set}_${name}_alignment.sam > \${k2}_${name}.bow1_\${rna_set}  2>&1
                    mv unmapped/${name}.unmapped_1.fastq unmapped/${name}.unmapped.1.fastq
                    mv unmapped/${name}.unmapped_2.fastq unmapped/${name}.unmapped.2.fastq
                fi
            else
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un  unmapped/${name}.unmapped.fastq -U ${name}.fastq --al ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam > \${k2}_${name}.bow_\${rna_set}  2>&1
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}  \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq  ${name}.fastq  -S \${rna_set}_${name}_alignment.sam > \${k2}_${name}.bow1_\${rna_set}  2>&1
                    
                fi
            fi
            echo "INFO: samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam"
            samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam
            if [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_tmp0.bam
                echo "INFO: samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam"
                samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam  # Remove unmapped reads
                if [ "${mate}" == "pair" ]; then
                    echo "# unique mapped reads: \$(samtools view -f 0x40 -F 0x4 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                else
                    echo "# unique mapped reads: \$(samtools view -F 0x40 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                fi
            fi
            if [ "${mate}" == "pair" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_alignment.tmp1.bam
                echo "INFO: samtools sort -n -o \${rna_set}_${name}_alignment.tmp2 \${rna_set}_${name}_alignment.tmp1.bam"
                samtools sort -n -o \${rna_set}_${name}_alignment.tmp2.bam \${rna_set}_${name}_alignment.tmp1.bam 
                echo "INFO: samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam"
                samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam
                rm \${rna_set}_${name}_alignment.tmp1.bam \${rna_set}_${name}_alignment.tmp2.bam
            fi
            echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam"
            samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam 
            echo "INFO: samtools index \${rna_set}@${name}_sorted.bam"
            samtools index \${rna_set}@${name}_sorted.bam
            
            if [ "${remove_duplicates}" == "yes" ]; then
                ## check read header whether they have UMI tags which are separated with underscore.(eg. NS5HGY:2:11_GTATAACCTT)
                umiCheck=\$(samtools view \${rna_set}@${name}_sorted.bam |head -n 1 | awk 'BEGIN {FS="\\t"}; {print \$1}' | awk 'BEGIN {FS=":"}; \$NF ~ /_/ {print \$NF}')
                
                # based on remove_duplicates_based_on_UMI_after_mapping
                if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes" -a ! -z "\$umiCheck" ]; then
                    echo "INFO: umi_mark_duplicates.py will be executed for removing duplicates from bam file"
                    echo "python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4"
                    python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4
                else
                    echo "INFO: Picard MarkDuplicates will be executed for removing duplicates from bam file"
                    if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes"  ]; then
                        echo "WARNING: Read header have no UMI tags which are separated with underscore. Picard MarkDuplicates will be executed to remove duplicates from alignment file (bam) instead of remove_duplicates_based_on_UMI_after_mapping."
                    fi
                    echo "INFO: picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam"
                    picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam 
                fi
                #get duplicates stats (read the sam flags)
                samtools flagstat \${rna_set}@${name}_sorted.deumi.sorted.bam > \${k2}@\${rna_set}@${name}_duplicates_stats.log
                #remove alignments marked as duplicates
                samtools view -b -F 0x400 \${rna_set}@${name}_sorted.deumi.sorted.bam > \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup
                #sort deduplicated files by chrom pos
                echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup"
                samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup 
                samtools index \${rna_set}@${name}_sorted.dedup.bam
                #get flagstat after dedup
                echo "##After Deduplication##" >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                samtools flagstat \${rna_set}@${name}_sorted.dedup.bam >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
            fi
            
        
            for file in unmapped/*; do mv \$file \${file/.unmapped/}; done ##remove .unmapped from filename
            if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                grep -v Warning \${k2}_${name}.bow_\${rna_set} > ${name}.tmp
                mv ${name}.tmp \${k2}_${name}.bow_\${rna_set}
                cp \${k2}_${name}.bow_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                cp \${k2}_${name}.bow1_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                cp \${k2}_${name}.star_\${rna_set} ./../bowfiles/.
            fi
            cd ..
        else
            echo "WARNING: \${indexesListAr[\$k-1]} Mapping skipped. File not found."
            cd unmapped 
            ln -s \${wrkDir}/\${rna_set}/*fastq .
            cd ..
            cd ..
        fi
    done
    cd final_reads && ln -s \${wrkDir}/\${prev}/* .
else 
    mv ${reads} final_reads/.
fi
"""
}
}


params.run_Split_Fastq =  "no"  //* @dropdown @options:"yes","no" @show_settings:"SplitFastq"
readsPerFile = params.SplitFastq.readsPerFile
//Since splitFastq operator requires flat file structure, first convert grouped structure to flat, execute splitFastq, and then return back to original grouped structure
//.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

//Mapping grouped read structure to flat structure
flatPairsClosure = {row -> if(row[1] instanceof Collection) {
        if (row[1][1]){
            tuple(row[0], file(row[1][0]), file(row[1][1]))
        } else {
            tuple(row[0], file(row[1][0]))
        }
    } else {
        tuple(row[0], file(row[1]))
    }
}

//Mapping flat read structure to grouped read structure
groupPairsClosure = {row -> tuple(row[0], (row[2]) ? [file(row[1]), file(row[2])] : [file(row[1])])}

// if mate of split process different than rest of the pipeline, use "mate_split" as input parameter. Otherwise use default "mate" as input parameter
mateParamName = (params.mate_split) ? "mate_split" : "mate"
splitFastqParams = ""
if (params[mateParamName] != "pair"){
    splitFastqParams = [by: readsPerFile, file:true]
}else {
    splitFastqParams = [by: readsPerFile, pe:true, file:true]
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!(params.run_Split_Fastq == "yes")){
g124_32_reads_g_68.into{g_68_reads_g127_13}
} else {


process SplitFastq {

input:
 val mate from g_122_mate_g_68
 set val(name), file(reads) from g124_32_reads_g_68.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

output:
 set val(name), file("split/*")  into g_68_reads_g127_13

when:
params.run_Split_Fastq == "yes"

script:
"""    
mkdir -p split
mv ${reads} split/.
"""
}
}


g125_6_genomeIndexPath_g127_13= g125_6_genomeIndexPath_g127_13.ifEmpty(file('Bowtie2Index', type: 'any')) 

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 18
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1500
    $CPU  = 3
    $MEMORY = 18
    $QUEUE = "long"
}
//* platform
//* autofill

process Bowtie2_Module_Map_Bowtie2 {

input:
 set val(name), file(reads) from g_68_reads_g127_13
 val mate from g_122_mate_g127_13
 val Bowtie2Index from g125_6_genomeIndexPath_g127_13

output:
 set val(name), file("${newName}.bow")  into g127_13_bowfiles_g127_10, g127_13_bowfiles_g_79
 set val(name), file("${newName}_alignment.bam")  into g127_13_bam_file_g127_15

when:
(params.run_Bowtie2 && (params.run_Bowtie2 == "yes")) || !params.run_Bowtie2

script:
Map_Bowtie2_parameters = params.Bowtie2_Module_Map_Bowtie2.Map_Bowtie2_parameters

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = "";

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

""" 
    if [ "${mate}" == "pair" ]; then
        bowtie2 -x ${params.bowtie2_index} ${Map_Bowtie2_parameters} --no-unal  -1 ${file1} -2 ${file2} -S ${newName}_alignment.sam > ${newName}.bow 2>&1
    else
        bowtie2 -x ${params.bowtie2_index} ${Map_Bowtie2_parameters}  -U ${file1} -S ${newName}_alignment.sam > ${newName}.bow 2>&1
    fi
    grep -v Warning ${newName}.bow > ${newName}.tmp
    mv  ${newName}.tmp ${newName}.bow 
    samtools view -bS ${newName}_alignment.sam > ${newName}_alignment.bam 
    rm -f ${newName}_alignment.sam
"""


}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Bowtie2_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bam$/) "bowtie2/$filename"
}

input:
 set val(oldname), file(bamfiles) from g127_13_bam_file_g127_15.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g127_15_merged_bams
 set val(oldname), file("*_sorted*bai")  into g127_15_bam_index
 set val(oldname), file("*_sorted*bam")  into g127_15_sorted_bam_g128_21

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -O bam -T !{oldname} -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -T !{oldname} -O bam -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

params.run_Remove_Multimappers_with_Samtools =  "no"  //* @dropdown @options:"yes","no" @show_settings:"Remove_Multimappers"

if (!(params.run_Remove_Multimappers_with_Samtools == "yes")){
g127_15_sorted_bam_g128_21.into{g128_21_mapped_reads_g128_22}
} else {


process ChIP_Module_Remove_Multimappers_with_Samtools {

input:
 set val(name), file(bam) from g127_15_sorted_bam_g128_21

output:
 set val(name), file("bam/${name}.bam")  into g128_21_mapped_reads_g128_22

when:
params.run_Remove_Multimappers_with_Samtools == "yes" 

script:
MAPQ_quality = params.ChIP_Module_Remove_Multimappers_with_Samtools.MAPQ_quality
"""
mkdir bam
samtools view -hb -q ${MAPQ_quality} ${bam} > ${name}_unique.bam
mv ${name}_unique.bam bam/${name}.bam
"""

}
}


mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Dedup_Bam_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "sequential_mapping_counts/$filename"
}

input:
 file bam from g124_32_bam_file_g124_27.collect()
 file index from g124_32_bam_index_g124_27.collect()

output:
 file "*.counts.tsv"  into g124_27_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        print $3;
        if ($3 eq ".dedup"){
            $dedup = "dedup.";
        }
        push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
        push @{$all_files{$mapper}}, $file;
}


open OUT, ">header.tsv";
print OUT join ("\\t", "id","len",@header),"\\n";
close OUT;

foreach my $key (sort keys %all_files) {  
    my @array = @{ $all_files{$key} };  
        unless (-e "$indexHash{$key}.bed") {
        print "2: bed not found run makeBed\\n";
            if (-e "$indexHash{$key}.fa") { 
                makeBed("$indexHash{$key}.fa", $key, "$indexHash{$key}.bed");
            } elsif(-e "$indexHash{$key}.fasta"){
                makeBed("$indexHash{$key}.fasta", $key, "$indexHash{$key}.bed");
            }
        }
    
        my $bamFiles = join ' ', @array;
        print "bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp\\n";
        `bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp`;
        my $iniResColumn = int(countColumn("$indexHash{$key}.bed")) + 1;
        `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key.${dedup}counts.tmp> $key.${dedup}counts.tsv`;
        `sort -k3,3nr $key.${dedup}counts.tsv>$key.${dedup}sorted.tsv`;
        `cat header.tsv $key.${dedup}sorted.tsv> $key.${dedup}counts.tsv`;
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''


}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Bam_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "sequential_mapping_counts/$filename"
}

input:
 file bam from g124_32_bam_file_g124_23.collect()
 file index from g124_32_bam_index_g124_23.collect()

output:
 file "*.counts.tsv"  into g124_23_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        print $3;
        if ($3 eq ".dedup"){
            $dedup = "dedup.";
        }
        push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
        push @{$all_files{$mapper}}, $file;
}


open OUT, ">header.tsv";
print OUT join ("\\t", "id","len",@header),"\\n";
close OUT;

foreach my $key (sort keys %all_files) {  
    my @array = @{ $all_files{$key} };  
        unless (-e "$indexHash{$key}.bed") {
        print "2: bed not found run makeBed\\n";
            if (-e "$indexHash{$key}.fa") { 
                makeBed("$indexHash{$key}.fa", $key, "$indexHash{$key}.bed");
            } elsif(-e "$indexHash{$key}.fasta"){
                makeBed("$indexHash{$key}.fasta", $key, "$indexHash{$key}.bed");
            }
        }
    
        my $bamFiles = join ' ', @array;
        print "bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp\\n";
        `bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp`;
        my $iniResColumn = int(countColumn("$indexHash{$key}.bed")) + 1;
        `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key.${dedup}counts.tmp> $key.${dedup}counts.tsv`;
        `sort -k3,3nr $key.${dedup}counts.tsv>$key.${dedup}sorted.tsv`;
        `cat header.tsv $key.${dedup}sorted.tsv> $key.${dedup}counts.tsv`;
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''


}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Deduplication_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /deduplication_summary.tsv$/) "sequential_mapping_summary/$filename"
}

input:
 file flagstat from g124_32_log_file_g124_30.collect()
 val mate from g_122_mate_g124_30

output:
 file "deduplication_summary.tsv"  into g124_30_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Sequential_Mapping_Module_Sequential_Mapping_Summary {

input:
 set val(name), file(bowfile) from g124_32_bowfiles_g124_26
 val mate from g_122_mate_g124_26
 val filtersList from g124_32_filter_g124_26

output:
 file '*.tsv'  into g124_26_outputFileTSV_g124_13
 val "sequential_mapping_sum"  into g124_26_name_g124_13

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', "!{name}.tsv");
print $fh "Sample\\tGroup\\tTotal Reads\\tReads After Sequential Mapping\\tUniquely Mapped\\tMultimapped\\tMapped\\n";
my @bowArray = split(' ', "!{bowfile}");
my $group= "\\t";
my @filterArray = (!{filtersList});
foreach my $bowitem(@bowArray) {
    # get mapping id
    my @bowAr = $bowitem.split("_");
    $bowCount = $bowAr[0] + -1;
    # if bowfiles ends with underscore (eg. bow_rRNA), parse rRNA as a group.
    my ($RDS_In, $RDS_After, $RDS_Uniq, $RDS_Multi, $ALGN_T, $a, $b, $aPer, $bPer)=(0, 0, 0, 0, 0, 0, 0, 0, 0);
    if ($bowitem =~ m/bow_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN, $bowitem);
        my $i = 0;
        while(my $line=<IN>){
            chomp($line);
            $line=~s/^ +//;
            my @arr=split(/ /, $line);
            $RDS_In=$arr[0] if ($i=~/^1$/);
            # Reads After Filtering column depends on filtering type
            if ($i == 2){
                if ($filterArray[$bowCount] eq "Yes"){
                    $RDS_After=$arr[0];
                } else {
                    $RDS_After=$RDS_In;
                }
            }
            if ($i == 3){
                $a=$arr[0];
                $aPer=$arr[1];
                $aPer=~ s/([()])//g;
                $RDS_Uniq=$arr[0];
            }
            if ($i == 4){
                $b=$arr[0];
                $bPer=$arr[1];
                $bPer=~ s/([()])//g;
                $RDS_Multi=$arr[0];
            }
            $ALGN_T=($a+$b);
            $i++;
        }
        close(IN);
    } elsif ($bowitem =~ m/star_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $bowitem | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $bowitem | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($multimapped);
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = $RDS_Uniq+$RDS_Multi;
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    } elsif ($bowitem =~ m/bow1_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		my $uniqAligned;
		chomp($inputCount = `cat $bowitem | grep '# reads processed:' | awk '{sum+=\\$4} END {print sum}'`);
		chomp($aligned = `cat $bowitem | grep '# reads with at least one reported alignment:' | awk '{sum+=\\$9} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep '# unique mapped reads:' | awk '{sum+=\\$5} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($aligned) -int($uniqAligned);
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = int($aligned);
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    }
    
    print $fh "!{name}\\t$group$RDS_In\\t$RDS_After\\t$RDS_Uniq\\t$RDS_Multi\\t$ALGN_T\\n";
}
close($fh);



'''

}


process Sequential_Mapping_Module_Merge_TSV_Files {

input:
 file tsv from g124_26_outputFileTSV_g124_13.collect()
 val outputFileName from g124_26_name_g124_13.collect()

output:
 file "${name}.tsv"  into g124_13_outputFileTSV_g124_14

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process Sequential_Mapping_Module_Sequential_Mapping_Short_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequential_mapping_short_sum.tsv$/) "sequential_mapping_summary/$filename"
	else if (filename =~ /sequential_mapping_detailed_sum.tsv$/) "sequential_mapping_summary/$filename"
}

input:
 file mainSum from g124_13_outputFileTSV_g124_14

output:
 file "sequential_mapping_short_sum.tsv"  into g124_14_outputFileTSV_g_115
 file "sequential_mapping_detailed_sum.tsv"  into g124_14_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols_short;
my @seen_cols_detailed;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @h) = ( split("\\t", $line1) );
        my $totalHeader = $h[1];
        my $afterFilteringHeader = $h[2];
        my $uniqueHeader = $h[3];
        my $multiHeader = $h[4];
        my $mappedHeader = $h[5];
        push(@seen_cols_short, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_short; #Total reads Header
        push(@seen_cols_detailed, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_detailed; #Total reads Header

        my $n=0;
        while (my $line=<IN>) {
                
                chomp($line);
                my ( $ID, @fields ) = ( split("\\t", $line) ); 
                #SHORT
                push(@seen_cols_short, $fields[0]) unless grep{$_ eq $fields[0]} @seen_cols_short; #mapped element header
                $all_rows{$ID}{$fields[0]} = $fields[5];#Mapped Reads
                #Grep first line $fields[1] as total reads.
                if (!exists $all_rows{$ID}{$totalHeader}){    
                        $all_rows{$ID}{$totalHeader} = $fields[1];
                } 
                $all_rows{$ID}{$afterFilteringHeader} = $fields[2]; #only use last entry
                #DETAILED
                $uniqueHeadEach = "$fields[0] (${uniqueHeader})";
                $multiHeadEach = "$fields[0] (${multiHeader})";
                $mappedHeadEach = "$fields[0] (${mappedHeader})";
                push(@seen_cols_detailed, $mappedHeadEach) unless grep{$_ eq $mappedHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $uniqueHeadEach) unless grep{$_ eq $uniqueHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $multiHeadEach) unless grep{$_ eq $multiHeadEach} @seen_cols_detailed;
                $all_rows{$ID}{$mappedHeadEach} = $fields[5];
                $all_rows{$ID}{$uniqueHeadEach} = $fields[3];
                $all_rows{$ID}{$multiHeadEach} = $fields[4];
    }
    close IN;
    push(@seen_cols_short, $afterFilteringHeader) unless grep{$_ eq $afterFilteringHeader} @seen_cols_short; #After filtering Header
}


#print Dumper \\%all_rows;
#print Dumper \\%seen_cols_short;

printFiles("sequential_mapping_short_sum.tsv",@seen_cols_short,);
printFiles("sequential_mapping_detailed_sum.tsv",@seen_cols_detailed);


sub printFiles {
    my($summary, @cols_to_print) = @_;
    
    open OUT, ">$summary";
    print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
    foreach my $key ( keys %all_rows ) { 
        print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
        }
        close OUT;
}

'''


}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Module_featureCounts_Prep {

input:

output:
 val run_params  into g126_125_run_parameters_g126_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_featureCounts_Prep.run_parameters

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
"""
"""

}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Bowtie2_Module_Bowtie_Summary {

input:
 set val(name), file(bowfile) from g127_13_bowfiles_g127_10.groupTuple()
 val mate from g_122_mate_g127_10

output:
 file '*.tsv'  into g127_10_outputFileTSV_g127_11
 val "bowtie_sum"  into g127_10_name_g127_11

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', "!{name}.tsv");
print $fh "Sample\\tTotal Reads\\tUnique Reads Aligned (Bowtie2)\\tMultimapped Reads Aligned (Bowtie2)\\n";
my @bowArray = split(' ', "!{bowfile}");
my ($RDS_T, $RDS_C1, $RDS_C2)=(0, 0, 0);
foreach my $bowitem(@bowArray) {
    # get mapping id
    open(IN, $bowitem);
    my $i = 0;
    while(my $line=<IN>)
    {
        chomp($line);
        $line=~s/^ +//;
        my @arr=split(/ /, $line);
        $RDS_T+=$arr[0] if ($i=~/^1$/);
        if ($i == 3){
            $RDS_C1+=$arr[0];
        }
        if ($i == 4){
            $RDS_C2+=$arr[0];
        }
        $i++;
    }
    close(IN);
}
print $fh "!{name}\\t$RDS_T\\t$RDS_C1\\t$RDS_C2\\n";
close($fh);



'''

}


process Bowtie2_Module_Merge_TSV_Files {

input:
 file tsv from g127_10_outputFileTSV_g127_11.collect()
 val outputFileName from g127_10_name_g127_11.collect()

output:
 file "${name}.tsv"  into g127_11_outputFileTSV_g_115

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Remove_Multimappers_with_Picard && (params.run_Remove_Multimappers_with_Picard == "yes")) || !params.run_Remove_Multimappers_with_Picard)){
g128_21_mapped_reads_g128_22.into{g128_22_mapped_reads_g128_25; g128_22_mapped_reads_g126_121; g128_22_mapped_reads_g126_122; g128_22_mapped_reads_g126_123; g128_22_mapped_reads_g126_124; g128_22_mapped_reads_g126_126}
g128_22_publish = Channel.empty()
g128_22_log_file_g128_23 = Channel.empty()
} else {


process ChIP_Module_Picard_MarkDuplicates {

input:
 set val(name), file(bam) from g128_21_mapped_reads_g128_22

output:
 set val(name), file("bam/${name}.bam")  into g128_22_mapped_reads_g128_25, g128_22_mapped_reads_g126_121, g128_22_mapped_reads_g126_122, g128_22_mapped_reads_g126_123, g128_22_mapped_reads_g126_124, g128_22_mapped_reads_g126_126
 set val(name), file("${name}*")  into g128_22_publish
 file "*_duplicates_stats.log"  into g128_22_log_file_g128_23

when:
(params.run_Remove_Multimappers_with_Picard && (params.run_Remove_Multimappers_with_Picard == "yes")) || !params.run_Remove_Multimappers_with_Picard     

script:
"""
mkdir bam
picard MarkDuplicates OUTPUT=${name}_dedup.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=${bam} > ${name}_picard.log 

#get duplicates stats (read the sam flags)
samtools flagstat ${name}_dedup.bam > ${name}@Reads@${name}_duplicates_stats.log
#remove alignments marked as duplicates
samtools view -b -F 0x400 ${name}_dedup.bam > ${name}_dedup.bam.x_dup
#sort deduplicated files by chrom pos
samtools sort -o ${name}_sorted.dedup.bam ${name}_dedup.bam.x_dup 
mv ${name}_sorted.dedup.bam bam/${name}.bam
#get properly paired log stats after dedup
echo "##After Deduplication##" >> ${name}@Reads@${name}_duplicates_stats.log
samtools flagstat bam/${name}.bam >> ${name}@Reads@${name}_duplicates_stats.log
"""
}
}


params.gtf =  ""  //* @input


process BAM_Analysis_Module_featureCounts {

input:
 set val(name), file(bam) from g128_22_mapped_reads_g126_126
 val paired from g_122_mate_g126_126
 each run_params from g126_125_run_parameters_g126_126

output:
 file "*"  into g126_126_outputFileTSV_g126_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Module_summary_featureCounts {

input:
 file featureCountsOut from g126_126_outputFileTSV_g126_117.collect()

output:
 file "*_featureCounts.tsv"  into g126_117_outputFile
 file "*_featureCounts.sum.tsv"  into g126_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
    }
}

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

params.bed =  ""  //* @input

process BAM_Analysis_Module_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc/$filename"
}

input:
 set val(name), file(bam) from g128_22_mapped_reads_g126_122

output:
 file "result/*.out"  into g126_122_outputFileOut_g126_95, g126_122_outputFileOut_g_79

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary/$filename"
}

input:
 file rseqcOut from g126_122_outputFileOut_g126_95.collect()
 val mate from g_122_mate_g126_95

output:
 file "*.tsv"  into g126_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Module_Picard {

input:
 set val(name), file(bam) from g128_22_mapped_reads_g126_121

output:
 file "*_metrics"  into g126_121_outputFileOut_g126_82
 file "results/*.pdf"  into g126_121_outputFilePdf_g126_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard/$filename"
}

input:
 file picardOut from g126_121_outputFileOut_g126_82.collect()
 val mate from g_122_mate_g126_82
 file picardPdf from g126_121_outputFilePdf_g126_82.collect()

output:
 file "*.tsv"  into g126_82_outputFileTSV
 file "results/*.pdf"  into g126_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}

params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_Module_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "ucsc_bigwig/$filename"
}

input:
 set val(name), file(bam) from g128_22_mapped_reads_g126_124

output:
 file "*.bw"  into g126_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = ''
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index "+ name+"_sorted.bam "
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${name}_sorted.bam -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_Module_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_IGV_BAM2TDF_converter.igv_window_size

params.genome =  ""  //* @input

process BAM_Analysis_Module_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igv_tdf/$filename"
}

input:
 val mate from g_122_mate_g126_123
 set val(name), file(bam) from g128_22_mapped_reads_g126_123

output:
 file "*.tdf"  into g126_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = ''
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index "+ name+"_sorted.bam "
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${name}_sorted.bam ${name}.tdf ${params.genome}
"""
}

macs2_callpeak_parameters = params.ChIP_Module_ChIP_Prep.macs2_callpeak_parameters
peak_calling_type = params.ChIP_Module_ChIP_Prep.peak_calling_type
band_width = params.ChIP_Module_ChIP_Prep.band_width
bedtoolsCoverage_Parameters = params.ChIP_Module_ChIP_Prep.bedtoolsCoverage_Parameters
compare_Custom_Bed = params.ChIP_Module_ChIP_Prep.compare_Custom_Bed
output_prefix = params.ChIP_Module_ChIP_Prep.output_prefix
sample_prefix = params.ChIP_Module_ChIP_Prep.sample_prefix
input_prefix = params.ChIP_Module_ChIP_Prep.input_prefix
//* @array:{output_prefix,sample_prefix,input_prefix} @multicolumn:{output_prefix,sample_prefix,input_prefix},{macs2_callpeak_parameters,peak_calling_type,band_width,bedtoolsCoverage_Parameters}
samplehash = [:]
inputhash = [:]
output_prefix.eachWithIndex { key, i -> inputhash[key] = input_prefix[i] }
output_prefix.eachWithIndex { key, i -> samplehash[key] = sample_prefix[i] }

process ChIP_Module_ChIP_Prep {

input:
 val mate from g_122_mate_g128_25
 set val(name), file(bam) from g128_22_mapped_reads_g128_25

output:
 file "bam/*.bam"  into g128_25_bam_file_g128_9
 val output_prefix  into g128_25_name_g128_9

when:
(params.run_ChIP_MACS2 && (params.run_ChIP_MACS2 == "yes")) || !params.run_ChIP_MACS2

script:
"""
mkdir -p bam
mv ${bam} bam/${name}.bam
"""
}


process ChIP_Module_ChIP_MACS2 {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /bam\/.*.bam$/) "chip/$filename"
	else if (filename =~ /${name}.*$/) "chip/$filename"
}

input:
 val mate from g_122_mate_g128_9
 file bam from g128_25_bam_file_g128_9.collect()
 val name from g128_25_name_g128_9.unique().flatten()

output:
 val compare_bed  into g128_9_compare_bed_g128_27
 file "*${peak_calling_type}Peak"  into g128_9_bed_g128_10, g128_9_bed_g132_1
 set val(name), file("bam/*.bam")  into g128_9_bam_file_g128_10, g128_9_bam_file_g128_27
 file "${name}*"  into g128_9_resultsdir_g_79
 val name  into g128_9_name

script:
genomeSizeText = ""
if (params.genome_build.contains("mouse")){
    genomeSizeText = "-g mm"
} else if (params.genome_build.contains("human")){
    genomeSizeText = "-g hs"
}

if (peak_calling_type == "narrow"){
    peakcallingType = ""
} else if (peak_calling_type == "broad"){
    peakcallingType = "--broad"
}

compare_bed = "merged.bed"
compare_Custom_Bed = compare_Custom_Bed.trim();
if (compare_Custom_Bed != ""){
    compare_bed = compare_Custom_Bed
}
inputsList = inputhash[name] 
samplesList = samplehash[name]

"""
echo ${samplesList}
echo ${inputsList}
echo $name
mkdir -p bam

#samplesList
samplesList="\$(echo -e "${samplesList}" | tr -d '[:space:]')" 
IFS=',' read -ra eachSampleAr <<< "\${samplesList}"
numSamples=\${#eachSampleAr[@]}
eachSampleArBam=( "\${eachSampleAr[@]/%/.bam }" )
sample_set=\${eachSampleArBam[@]}
bam_set=\${eachSampleArBam[@]}

#inputsList
input_set=""
inputsList="\$(echo -e "${inputsList}" | tr -d '[:space:]')" 
if [ "\${inputsList}" != "" ]; then
    IFS=',' read -ra eachInputAr <<< "\${inputsList}"
    eachInputArbam=( "\${eachInputAr[@]/%/.bam }" )
    input_set="-c \${eachInputArbam[@]}" 
    
fi
echo \${eachSampleArBam[@]}

macs2 callpeak --bw ${band_width} -t \${sample_set} \${input_set} -n ${name} ${genomeSizeText} ${macs2_callpeak_parameters} ${peakcallingType}

#bam files
if [ "\$numSamples" -gt "1" ]; then
    samtools merge bam/${name}.bam \$bam_set
else 
    rsync -a  \$bam_set bam/${name}.bam
fi

"""
}

params.run_Scripture =  "no"  //* @dropdown @options:"yes","no" @show_settings:"Scripture_peakrescore"
params.peakrescore_path =  ""  //* @input
params.peakrescore_class_path =  ""  //* @input

if (!(params.run_Scripture == "yes")){
g128_9_bed_g128_10.into{g128_10_bed_g128_26}
} else {


process ChIP_Module_Scripture_peakrescore {

input:
 file bed from g128_9_bed_g128_10
 set val(name), file(bam) from g128_9_bam_file_g128_10

output:
 file "${name}_trim.bed"  into g128_10_bed_g128_26

when:
params.run_Scripture == "yes"

script:
window = params.ChIP_Module_Scripture_peakrescore.window
trimFraction = params.ChIP_Module_Scripture_peakrescore.trimFraction
windowText = (window.toString() != "") ? "-window ${window}" : ""
trimFractionText = (trimFraction.toString() != "") ? "-trimFraction ${trimFraction}" : ""
"""
samtools index ${bam}
cat ${bed} | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5}' > ${name}_clean 
java -cp ${params.peakrescore_path}:${params.peakrescore_class_path} peaks.PeakTrim -task trimByFractionOfScale -in ${name}_clean -libAlignment ${bam}  $windowText $trimFractionText -out ${name}_trim.bed 
"""

}
}



process ChIP_Module_bed_merge {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /merged.bed$/) "chip/$filename"
}

input:
 file bed from g128_10_bed_g128_26.collect()

output:
 file "merged.bed"  into g128_26_bed_g128_27

"""
 cat ${bed} | cut -f -6 | bedtools sort -i stdin | bedtools slop -i stdin -b 100 -g ${params.genome_sizes} | bedtools merge -i stdin | awk '{print \$0"\t"\$1"_"\$2"_"\$3}' > merged.bed

"""
}

params.homer_dir =  ""  //* @input

process Homer_Find_Motif_Module_homer_download_install {

input:
 val run from g_119_run_Homer_g132_7

output:
 val run  into g132_7_run_process_g132_1

when:
run == "yes"

script:
"""
export PATH=\$PATH:${params.homer_dir}/bin/
if [ ! -d "${params.homer_dir}/data" ]; then
  echo "${params.homer_dir}/data not found"
  wget http://homer.ucsd.edu/homer/configureHomer.pl
  mkdir -p ${params.homer_dir}
  cp configureHomer.pl ${params.homer_dir}/.
  perl ${params.homer_dir}/configureHomer.pl -install
fi
if [ ! -d "${params.homer_dir}/data/genomes/${params.build}" ]; then
  echo "${params.build} will be installed"
  perl ${params.homer_dir}/configureHomer.pl -install ${params.build}
fi

"""

}

params.homer_dir =  ""  //* @input

process Homer_Find_Motif_Module_homer_find_Motifs_Genome {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /homer_${name}$/) "motif_finder_result_on_chip_peaks/$filename"
}

input:
 file bed from g128_9_bed_g132_1
 val run from g132_7_run_process_g132_1

output:
 file "homer_${name}"  into g132_1_resultsdir_g132_5

when:
run == "yes"

script:
name = bed.baseName
size = params.Homer_Find_Motif_Module_homer_find_Motifs_Genome.size
"""
mkdir preparsed
export PATH=\$PATH:${params.homer_dir}/bin/
findMotifsGenome.pl $bed ${_build} homer_${name} -size $size -preparsedDir preparsed
"""


}


process Homer_Find_Motif_Module_homer_summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /homerReport.*$/) "motif_finder_summary_on_chip_peaks/$filename"
}

input:
 file dir from g132_1_resultsdir_g132_5.collect()

output:
 file "homerReport*"  into g132_5_resultsdir

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $indir = $ENV{'PWD'};

opendir D, $indir or die "Could not open $indir";
my @alndirs = sort { $a cmp $b } grep /^homer_/, readdir(D);
closedir D;

`mkdir -p homerReport`;

foreach my $d (@alndirs){
    my $dir = "${indir}/$d";
    my $libname=$d;
    $libname=~s/homer_//;
    print "replace links in the html report with hidden and renamed folders";
    `cat $dir/homerResults.html | sed -e 's/.svg"/.html"/g' | sed -e 's:HREF="homerResults/:HREF=".${libname}homerResults/:g' | sed -e 's:<A HREF="knownResults.html">Known Motif Enrichment Results</A><BR/>::g' | sed -e 's:<A HREF="geneOntology.html">Gene Ontology Enrichment Results</A><BR/>::g' >> homerReport/homer_sum.html` if -e "$dir/homerResults.html";
    `cat $dir/knownResults.html | sed -e 's/.svg"/.html"/g' | sed -e 's:HREF="knownResults/:HREF=".${libname}knownResults/:g' | sed -e 's:HREF="knownResults.txt":HREF=".${libname}knownResults/knownResults.txt":g' | sed -e 's:<A HREF="homerResults.html">Homer <i>de novo</i> Motif Results</A><BR/>::g' | sed -e 's:<A HREF="geneOntology.html">Gene Ontology Enrichment Results</A><BR/>::g' >> homerReport/homer_sum.html` if -e "$dir/knownResults.html";
    `cp -r $dir/homerResults homerReport/.${libname}homerResults` if -d "$dir/homerResults";
    `cp -r $dir/knownResults homerReport/.${libname}knownResults` if -d "$dir/knownResults";
    `cp -r $dir/knownResults.txt homerReport/.${libname}knownResults/knownResults.txt` if -e "$dir/knownResults.txt";
    `find homerReport -name "*.html" -exec sed -i 's/.svg"/.html"/g' {} +`;
    `for fx in homerReport/.${libname}knownResults/*.svg; do mv \\$fx \\${fx\\%.svg}.html; done`;
    `for fx in homerReport/.${libname}homerResults/*.svg; do mv \\$fx \\${fx\\%.svg}.html; done`;
    `zip -r homerReport.zip homerReport`;
}


'''
}

params.run_FastQC =  "no"  //* @dropdown @options:"yes","no"



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(html|zip)$/) "fastqc/$filename"
}

input:
 val mate from g_122_mate_g123_3
 set val(name), file(reads) from g_1_reads_g123_3

output:
 file '*.{html,zip}'  into g123_3_FastQCout_g_79

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process MultiQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /multiqc_report.html$/) "multiQC/$filename"
}

input:
 file "fastqc/*" from g123_3_FastQCout_g_79.flatten().toList()
 file "sequential_mapping/*" from g124_32_bowfiles_g_79.flatten().toList()
 file "rseqc_hisat/*" from g126_122_outputFileOut_g_79.flatten().toList()
 file "macs/*" from g128_9_resultsdir_g_79.flatten().toList()
 file "bowtie/*" from g127_13_bowfiles_g_79.flatten().toList()

output:
 file "multiqc_report.html" optional true  into g_79_htmlout

"""
multiqc -e general_stats -d -dd 2 .
"""
}



process ChIP_Module_bedtools_coverage {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.sum.txt$/) "chip/$filename"
}

input:
 val compare_bed from g128_9_compare_bed_g128_27
 file bed from g128_26_bed_g128_27
 set val(name), file(bam) from g128_9_bam_file_g128_27

output:
 file "*.sum.txt"  into g128_27_outputFileTxt_g128_13

script:
bedtoolsCoverage_Parameters = params.ChIP_Module_bedtools_coverage.bedtoolsCoverage_Parameters
bedtoolsIntersect_Parameters = params.ChIP_Module_bedtools_coverage.bedtoolsIntersect_Parameters
"""
echo ${compare_bed}
if [ -s "${compare_bed}" ]; then 
    echo " bed file exists and is not empty "
        samtools view -H ${name}.bam | grep -P "@SQ\\tSN:" | sed 's/@SQ\\tSN://' | sed 's/\\tLN:/\\t/' > ${name}_chroms
        bedtools intersect -abam ${name}.bam -b ${compare_bed} > temp_${name}.bam
        bedtools sort -faidx ${name}_chroms -i ${compare_bed}  | bedtools coverage ${bedtoolsCoverage_Parameters} -a stdin -b temp_${name}.bam  > temp_${name}.bed
        # 'The number of features in B that overlapped the A interval' multiplied by 'fraction of bases in A that had non-zero coverage from features in B'.
        awk '{\$NF=\$(NF-3)*\$NF;print }' OFS="\\t" temp_${name}.bed | grep -v all > temp_${name}_hist.bed
        l=`awk '{print NF}' temp_${name}_hist.bed | head -1 | awk '{print \$1-4}'`
        k=`awk '{print NF}' temp_${name}_hist.bed | head -1`
        bedtools groupby -i temp_${name}_hist.bed -g 1-\$l -c \$k -o sum > ${name}.sum.txt
        #rm -rf temp_*

else
  echo " bed file does not exist, or is empty "
  touch ${name}_empty.sum.txt
fi
"""

}


process ChIP_Module_ATAC_CHIP_summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "chip_summary/$filename"
}

input:
 file file from g128_27_outputFileTxt_g128_13.collect()

output:
 file "*.tsv"  into g128_13_outputFile

shell:
'''
#!/usr/bin/env perl

my $indir = $ENV{'PWD'};

opendir D, $indir or die "Could not open $indir\n";
my @alndirs = sort { $a cmp $b } grep /.txt/, readdir(D);
closedir D;
    
my @a=();
my %b=();
my %c=();
my $i=0;
foreach my $d (@alndirs){ 
    my $file = "${indir}/$d";
    print $d."\n";
    my $libname=$d;
    $libname=~s/\\.sum\\.txt//;
    print $libname."\n";
    $i++;
    $a[$i]=$libname;
    open IN,"${indir}/$d";
    $_=<IN>;
    while(<IN>)
    {
        my @v=split; 
        $b{$v[3]}{$i}=$v[4];
    }
    close IN;
}
my $outfile="${indir}/"."sum_counts.tsv";
open OUT, ">$outfile";
print OUT "Feature";

for(my $j=1;$j<=$i;$j++) {
    print OUT "\t$a[$j]";
}
print OUT "\n";
    
foreach my $key (keys %b){
    print OUT "$key";
    for(my $j=1;$j<=$i;$j++){
        print OUT "\t$b{$key}{$j}";
    }
    print OUT "\n";
}
close OUT;
'''
}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process ChIP_Module_Picard_Deduplication_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /deduplication_summary.tsv$/) "picard_deduplication/$filename"
}

input:
 file flagstat from g128_22_log_file_g128_23.collect()
 val mate from g_122_mate_g128_23

output:
 file "deduplication_summary.tsv"  into g128_23_outputFileTSV_g_115

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}

g124_14_outputFileTSV_g_115= g124_14_outputFileTSV_g_115.ifEmpty(file('sequentialSum', type: 'any')) 
g127_11_outputFileTSV_g_115= g127_11_outputFileTSV_g_115.ifEmpty(file('hisatSum', type: 'any')) 
g123_11_outputFileTSV_g_115= g123_11_outputFileTSV_g_115.ifEmpty(file('adapterSum', type: 'any')) 
g123_21_outputFileTSV_g_115= g123_21_outputFileTSV_g_115.ifEmpty(file('trimmerSum', type: 'any')) 
g123_16_outputFileTSV_g_115= g123_16_outputFileTSV_g_115.ifEmpty(file('qualitySum', type: 'any')) 
g128_23_outputFileTSV_g_115= g128_23_outputFileTSV_g_115.ifEmpty(file('umiSum', type: 'any')) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /overall_summary.tsv$/) "summary/$filename"
}

input:
 file sequentialSum from g124_14_outputFileTSV_g_115
 file hisatSum from g127_11_outputFileTSV_g_115
 file adapterSum from g123_11_outputFileTSV_g_115
 file trimmerSum from g123_21_outputFileTSV_g_115
 file qualitySum from g123_16_outputFileTSV_g_115
 file umiSum from g128_23_outputFileTSV_g_115

output:
 file "overall_summary.tsv"  into g_115_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
my @order = ("adapter_removal","trimmer","quality","extractUMI","sequential_mapping","bowtie","star","hisat2","tophat2", "dedup","rsem");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
