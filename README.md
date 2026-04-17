# QRNG_With_QC
Work done as an Undergraduate Thesis in Physics by Lucas Hernández Bellón.

This repository contains all code used to extract, refine and study bitstrings of QRNG computed by Quantum Computers of IBM and AWS Cloud. 
Using a simple Hadamard-Transform Circuit, one can generate strings of n bits which ideally is random. The proposition of the TFG was to 
investigate to which point we can use the actual available QC to generate such randomness in the premise of *impredictibility*, *independance* and *uniformness* of bits in both the sequences and the space 
containing all sequences of length n {0,1}^n.

The repository is structured as follows:
  1. On *Root* there are the scripts, each focused in a specific task. We will talk about them later
  2. *Data_Processing* contains the study logs of the sequences obtained. The sequences per se are big enought to not being able to post here (I may link them with a dropbox account, I need to think about it)
  3. *Outputs* would contain any sequence obtained by the scripts by default. It can be changed.
  4. *NIST* contains all the tests of the *NIST800-22-1a Test Suite*, ported to Python by stevenang et.al.( https://github.com/stevenang/randomness_testsuite; On MIT License), but with some changes to focus only in the test's methods.   

## QuickStart Guide ##

For an easy Quickstart you need the following preliminaries:

1. An IBMQ Account and all configured to connect remotely to IBM Quantum Computers via token. Insert all needed in *Call_Computers.py* file, on *save_account* method of line 53.
   https://quantum.cloud.ibm.com/docs/en/guides/hello-world
3. An AWS Account. Configure the creedentials in your terminal as AWS guide indicates to be able to connect to their and their partner's computers.
   https://aws.amazon.com/es/blogs/quantum-computing/setting-up-your-local-development-environment-in-amazon-braket/
### Get Bitstrings ###
After that you can just run *main.py* on default variables and you will obtain some sequences as txts in the *Outputs* folder. Though we recommend to check first the number of qbits used, the length of the bitstring desired as well as the backends used.
If you dont need the AWS computers you can just switch off their methods in *main* and use only IBM Computers. 

### Analyse the Bitstrings ###
Given that you have some bitstrings stored as *main.py* produces them, you can use *main_Data.py* to analyse them and put the log results for each computer in the *Data_Processing* folder. You only need to execute *Data_Processing* given the strings' path as a variable. We left the results of our analysis in a subfolder called *TFG Results*, so feel free to check it before doing anything to see what to expect. 
