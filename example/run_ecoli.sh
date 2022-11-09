# check arguments
if [[ "$*" == *"--full"* ]]
then
    input=nearperfect-ecoli.10X.fa
else
    input=nearperfect-ecoli.100.fa
fi

# check availability of paftools.js
if [ -x "$(command -v paftools.js)" ]
then
    paftools=paftools.js
else
    if [ ! -f paftools.js ]
    then
       wget https://raw.githubusercontent.com/lh3/minimap2/master/misc/paftools.js
       chmod +x  paftools.js
    fi
    paftools=./paftools.js
fi


echo "hifimap -------------"


/usr/bin/time cargo run --release -- $input --reference ecoli.genome.fa --debug -k 8 -d 0.01 -l 16 -p mapquik -g 100 --threads 11
$paftools mapeval mapquik.paf

echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -ax asm20 -t 8 ecoli.genome.fa $input > minimap2.paf
$paftools mapeval minimap2.paf
