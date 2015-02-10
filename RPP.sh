#! /bin/bash

# written by Sishuo Wang from The University of British Columbia
# Please write e-mails to sishuowang@hotmail.ca if you have any question and/or suggestion. Your help is highly appreciated.
###################################################################

mnt3_sswang="/mnt/bay3/sswang"
fastq_dump="$mnt3_sswang/software/NGS/basic_processing/sratoolkit.2.4.2/bin/fastq-dump"
cutadapt_ruby="$mnt3_sswang/tools_program/NGS_scripts/basic/cutadapt.rb"
fastq_line3_del="$mnt3_sswang/tools_program/NGS_scripts/basic/fastq_line3_del.py"
fastq_detect="$mnt3_sswang/software/NGS/basic_processing/mini_tools/fastq_detect.pl"
#check_paired_or_single="$mnt3_sswang/tools_program/NGS_scripts/basic/check_paired_or_single_end.pl"

tophat2=tophat2
bowtie2=bowtie2
export PATH=$PATH:$mnt3_sswang/software/NGS/reads_map/bismark_package/bismark_v0.8.3/


###################################################################
function parse_sras(){
	for i in ${sras[@]}; do
		#$fastq_dump -O $SRA_examine_dir $i&
		unset a
		for j in `echo $(get_basename_corename $i)`; do
			a=(${a[@]} $j)
		done
		basename=${a[0]}
		corename=${a[1]}
		#tmp_fastq=$SRA_examine_dir/$corename".fastq"
		#proc_id=$!
		#sleep 1
		#kill -9 $proc_id
		#check_paired_or_single $tmp_fastq
		#paired_infos=(${paied_info[@]} $is_paired)
		corenames=(${corenames[@]} $corename) 
	done
}


function run_cutadapt(){
	unset cutadapt_fastqs
	declare -a cutadapt_fastqs
	if grep "[a-z]" <<< $cutadapt_ruby_args_content > /dev/null ; then
		cutadapt_ruby_args="$cutadapt_ruby_args_content"
	else
		cutadapt_ruby_args="no_args"
	fi
	for fastq in $@; do
		basename=`basename $fastq`
		output_basename=${basename/.fastq/.cutadapt.fastq}
		output=$cutadapt_outdir/$output_basename
		ruby $cutadapt_ruby -i $fastq --args $cutadapt_ruby_args -o $output --force
		cutadapt_fastqs=(${cutadapt_fastqs} $output)
	done
	echo ${cutadapt_fastqs[*]}
}


function map(){
	for index in ${!corenames[@]}; do
		echo "dump-split sra ......"
		corename=${corenames[$index]}
		sra=${sras[$index]}
		paired_info=${paired_infos[$index]}
		echo $corename

		outdir=$mapping_outdir/$corename
		mkdir $outdir

		$fastq_dump --split-3 -O $fastq_outdir $sra
		if [ -e $fastq_outdir/${corename}.fastq ]; then
			paired_info=0
		else
			paired_info=1
		fi

		if [ $paired_info == 1 ]; then
			fastq1=$fastq_outdir/"${corename}_1.fastq"
			fastq2=$fastq_outdir/"${corename}_2.fastq"
			fastqs=($fastq1 $fastq2)
			is_fastq_quality_phred64 $fastq1
			run_fastq_line3_del
			# run cutadapt
			if [ ! -z $is_cutadapt ]; then
				for i in $(run_cutadapt ${fastqs[@]}); do
					cutadapt_fastqs=(${cutadapt_fastqs[@]} $i)
				done
			fi

			fastq1=${cutadapt_fastqs[0]}
			fastq2=${cutadapt_fastqs[1]}
			case $mapper in
				tophat2)
					$tophat2 -p $cpu -o $outdir $genome_index $fastq1 $fastq2
					;;
				bowtie2)
					$bowtie2 -p $cpu -x $genome_index -1 $fastq1 -2 $fastq2 > $outdir/$corename.sam
					;;
				bismark)
					prepare_4_bismark_params
					bismark $non_directional_param --bowtie2 $phred64_param -o $outdir $bismark_genome_indir -1 $fastq1 -2 $fastq2
					fastq_basename=`basename $fastq1`
					cd $outdir
					ls
					bismark_methylation_extractor -s --comprehensive ${fastq_basename}_bismark_bt2_pe.sam
					cd -
					;;
			esac

		else
			fastq=$fastq_outdir/"$corename.fastq"
			fastqs=($fastq)
			is_fastq_quality_phred64 $fastq
			run_fastq_line3_del
			# run cutadapt
			if [ ! -z $is_cutadapt ]; then
				for i in $(run_cutadapt $fastq); do
					fastq=$i
				done
			fi
			case $mapper in
				tophat2)
					$tophat2 -p $cpu -o $outdir $genome_index $fastq
					;;
				bowtie2)
					$bowtie2 -p $cpu -x $genome_index -U $fastq > $outdir/$corename.sam
					;;
				bismark)
					prepare_4_bismark_params
					bismark $non_directional_param $phred64_param --bowtie2 -o $outdir $bismark_genome_indir $fastq
					fastq_basename=`basename $fastq`
					cd $outdir
					bismark_methylation_extractor -s --comprehensive ${fastq_basename}_bismark_bt2.sam
					cd -
					;;
			esac
		fi
	done
}


function is_fastq_quality_phred64(){
	local fastq_file
	fastq_file=$1
	quality_format="Phred33"
	line=`perl $fastq_detect $fastq_file 1000`
	if grep -P '  Illumina 1\.(3|5|13)\+\s+:  x' <<< $line | grep -P '  Solexa[ ]+:  x' <<< $line; then
		quality_format="Phred64"
	fi
	if [ $quality_format == "Phred64" ]; then
		is_phred64=true
	fi
}


function check_paired_or_single(){
	perl $check_paired_or_single $1
	if [ `perl $check_paired_or_single $1` -eq 1 ]; then
		is_paired=1
	else
		is_paired=0
	fi
	return $is_paired
}


function get_basename_corename(){
	basename=`basename $1`
	corename=${basename%%.sra}
	echo $basename $corename
}


function prepare_4_bismark_params(){
	if [ -z $is_strand_specific ]; then
		non_directional_param='--non_directional'
	fi
	if [ ! -z $is_phred64 ]; then
		phred64_param="--phred64-quals"
	fi
}


function run_fastq_line3_del(){
	for i in ${fastqs[@]}; do
		python $fastq_line3_del $i ${i}2
		mv ${i}2 $i
	done
}


function show_help(){
	echo $1
	basename=`basename $0`
	echo "bash $basename [Options]"
	cat <<EOF
Mandantory arguments:
--sra or
--fastq
--outdir
Optional arguments:
EOF
	exit
}


###################################################################
while [ $# -gt 0 ]; do
	case $1 in 
		--sra)
			sras=(${sras[@]} $2)
			shift
			;;
		--sra_dir)
			sra_dir=$2
			for i in $sra_dir/*.sra; do
				sras=(${sras[@]} $i)
			done
			shift
			;;
		--fastq)
			fastqs=(${fastqs[@]} $2)
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--mapper)
			mapper=$2
			shift
			;;
		--cutadapt)
			is_cutadapt=1
			;;
		--cutadapt_args)
			cutadapt_ruby_args_content=$2
			shift
			;;
		--cpu)
			cpu=$2
			shift
			;;
		--genome)
			genome=$2
			shift
			;;
		--genome_index)
			genome_index=$2
			shift
			;;
		--genome_indir|genome_indir_4_bismark)
			genome_indir_4_bismark=$2
			shift
			;;
		--bismark_genome_indir)
			bismark_genome_indir=$2
			shift
			;;			
		--strand_specific)
			is_strand_specific=1
			;;
		--clear|--force)
			force=1
			;;
		*)
			echo "Unknown option $1" >&2
			show_help
			;;
	esac
	shift
done


if [ -z $outdir ]; then
	show_help "Error! outdir should be specified with '--outdir'"
fi

if [ -d $outdir ]; then
	if [ $force ]; then
		rm -rf $outdir
	else
		show_help "outdir $outdir has already existed!"
	fi
fi


fastq_outdir=$outdir/"fastq"
mapping_outdir=$outdir/"mapping"
cutadapt_outdir=$outdir/"cutadapt"
bowtie2_genome_index_indir=$outdir/"bowtie2_genome_index"
mkdir -p $outdir
mkdir -p $fastq_outdir
mkdir -p $mapping_outdir
[ ! -z $is_cutadapt ] && mkdir -p $cutadapt_outdir


if [ $mapper == "bowtie2" -o $mapper == "tophat2" ]; then
	if [ -z $genome_index ]; then
		if [ ! -z $genome ]; then
			mkdir -p $bowtie2_genome_index_indir
			genome_basename=$bowtie2_genome_index_indir/`basename $genome`
			genome_basename=${genome_basename%.*}
			bowtie2-build $genome $genome_basename
			genome_index=$genome_basename
		else
			show_help "Error! Genome has to be given if genome_index is not given!"		
		fi
	fi
fi

if [ $mapper == "bismark" ]; then
	if [ -z $bismark_genome_indir ]; then
		if [ -z $genome_indir_4_bismark ]; then
			show_help "Error! bismark_genome_indir or genome_indir_4_bismark must be given if bismark is the mapper!"
		else
			bismark_genome_preparation --bowtie2 $genome_indir_4_bismark
			bismark_genome_indir=$genome_indir_4_bismark/
		fi
	fi
fi


[ -z $cpu ] && cpu=1


###################################################################
parse_sras

map


