#!/bin/sh

lines=100
seed=0; # unused currently. Need to replace shuf
out_dir=phylo.out
compat_args=""
iqtree3_args=""

ref=/nfs/srpipe_references/references/SARS-CoV-2/default/all/fasta/MN908947.3.fa

QDIR=/nfs/users/nfs_j/jkb/work/quantum/phylo
PATH=$QDIR:$PATH

help() {
    echo "Usage:covid_pipeline [-n lines] [-s seed] [-o out_dir] [-c compat_args] [-i iqtree3_args] covid.csv.xz"
}

while true
do
    case "$1" in
        "-h")
            help
            exit 0
            ;;
        "-n")
            lines=$2
            shift 2
            continue
            ;;
        "-s")
            seed=$2
            shift 2
            continue
            ;;
        "-o")
            out_dir=$2
            shift 2
            continue
            ;;
        "-c")
            compat_args=$2
            shift 2
            continue
            ;;
        "-i")
            iqtree3_args=$2
            shift 2
            continue
            ;;
        *)
            break
            ;;
    esac
done

if [ $# -eq 0 ]
then
    help
    exit
fi

tsv=$1

#--- Take a random selection of the input data
# Random lines with a seed.  Our alternative "shuf" command
shuf() {
    perl -e '$count=shift(@ARGV);
             srand(shift(@ARGV));
             while (<>) {push(@lines, $_)}
             for ($i=0;$i<$count;) {
                 $l = int(rand $#lines+.999);
                 next unless defined($lines[$l]);
                 print $lines[$l];
                 undef $lines[$l];
                 $i++;
             }' $1 $2
}

echo ":- Extracting $lines lines to FASTA"
mkdir -p $out_dir
echo "xz -d < $tsv | shuf $lines $seed | awk... > $out_dir/unaligned.fa"
xz -d < $tsv | shuf $lines $seed | \
    awk '{printf(">%04d_%s_%s\t%s\n%s\n",NR,$3,$2,$2,$NF)}' > $out_dir/unaligned.fa
cd $out_dir

#--- Produce a true tree by parsing the lineage data
echo ":- Computing lineage tree"
lineage_fa2tree.pl unaligned.fa > lineage.nwk


#--- Align sequences
echo ":- Aligning sequences with minimap2"
align_seqs.pl $ref unaligned.fa > aligned.fa
#echo ":- Aligning sequences with MAFFT"
#align_seqs_mafft.pl unaligned.fa > aligned.fa


#--- Build tree with compat
time_fmt='Elapsed %e %E\nCPU     %U\nSystem  %S\nMax RSS %M KB'
echo ":- Running compat $compat_args"
/usr/bin/time -f "$time_fmt" sh -c "eval compat.sh $compat_args aligned.fa aligned.compat.nwk 2> compat.out"

echo ":- Evaluating compat tree"
iqtree3 -rf lineage.nwk aligned.compat.nwk >rf.out1
tail -1 aligned.compat.nwk.rfdist


# --- Build tree with iqtree3
echo ":- Running iqtree3 $iqtree3_args -s"
/usr/bin/time -f "$time_fmt" sh -c "eval iqtree3 --redo -s aligned.fa $iqtree3_args > iqtree3.out"

echo ":- Evaluating iqtree3 tree"
iqtree3 -rf lineage.nwk aligned.fa.treefile >rf.out2
tail -1 aligned.fa.treefile.rfdist

echo ":- Trees:"
echo "$out_dir/lineage.nwk            Lineage computed tree"
echo "$out_dir/aligned.compat.nwk     Maximum Compat tree"
echo "$out_dir/aligned.fa.treefile    Maximum Likelihood tree"
echo
echo "Use https://itol.embl.de/ page to try viewing trees"

# 10    226 296
# 10    231 297
# 100   223 297
# 100   221 296
