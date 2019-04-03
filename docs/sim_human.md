# Generating simulated reads for Human genome

### Building PBSIM
First, We build PBSIM as follows:
```
cd codes/pbsim
./configure
make
cd ../..
```
### Obtaining sample real dataset
Next, we download a sample real dataset that will be used by PBSIM to generate simulated reads. We download the fastq files corresponding to a cell:
```
mkdir -p data/simulated
cd data/simulated
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.1.subreads.fastq
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.2.subreads.fastq
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.3.subreads.fastq
```

Then we merge downloaded files to get a single fastq file:
```
cat m130929_024849_42213_c100518541910000001823079209281311_s1_p0.{1..3}.subreads.fastq > m130929_024849_42213_c100518541910000001823079209281311_s1_p0.all.subreads.fastq
```

### Simulating reads
Then run pbsim using `run_pbsim.sh` script located in codes folder:
```
../../codes/run_pbsim.sh hg38.fa m130929_024849_42213_c100518541910000001823079209281311_s1_p0.all.subreads.fastq 1
```

Then we sample from the generated reads as follows:
```
mkdir ../selected
cd ../selected
../../codes/extract-pbsim/extract ../simulated/ 25000 1000 1000000 1> pb_human.fasta 2> pb_human.m5
```
In this case, the file `pb_human.fasta` contains 25000 reads at least 1000 bp long with proper PacBio style fasta header. In addition, the `pb_human.m5` contains base level alignments of each read to the reference genome in m5 format (check [here](https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format) for a description about m5 format).
