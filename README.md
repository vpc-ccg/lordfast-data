# lordFAST data
This repository contains:

* Codes, scripts, and instructions for generating simulated data for the lordFAST paper
* Code for evaluation of the alignments of the simulated reads

### NOTE: The codes are tested only on linux!
### NOTE: The submodule PBSIM is added from [https://github.com/yukiteruono/pbsim](https://github.com/yukiteruono/pbsim)

## Downloading the repository
```
git clone --recursive https://github.com/vpc-ccg/lordfast-extra.git
```

## Building PBSIM
```
cd codes/pbsim
./configure
make
```
Then the PBSIM binary is located at `codes/pbsim/src/pbsim`

## Generating simulated reads
For the **human genome** see [here](https://github.com/vpc-ccg/lordfast-extra/blob/master/docs/sim_human.md)
