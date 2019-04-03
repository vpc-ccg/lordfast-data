# Generating simulated reads for Human genome

The fastq file used for generating simulated reads is obtained as the following. First the fastq files corresponding to a cell were downloaded:
```
cd <rootDir>
mkdir -p data/human/simulated
cd data/human/simulated
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.1.subreads.fastq
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.2.subreads.fastq
wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.3.subreads.fastq
```

Then they were merged to get a single fastq file:
```
cat m130929_024849_42213_c100518541910000001823079209281311_s1_p0.1.subreads.fastq m130929_024849_42213_c100518541910000001823079209281311_s1_p0.2.subreads.fastq m130929_024849_42213_c100518541910000001823079209281311_s1_p0.3.subreads.fastq > m130929_024849_42213_c100518541910000001823079209281311_s1_p0.all.subreads.fastq
```

Then run pbsim using `run_pbsim.sh` script located in codes folder:
```
<rootDir>/codes/run_pbsim.sh hg38.fa m130929_024849_42213_c100518541910000001823079209281311_s1_p0.all.subreads.fastq 1
```

Then we sample from the generated reads and merge them:
```
mkdir ../selected
cd ../selected
<rootDir>/codes/extract-pbsim/extract ../simulated/ 25000 1000 1000000 1> pb_human.fasta 2> pb_human.m5
```
which generates 25000 reads with minimum length of 1000.

