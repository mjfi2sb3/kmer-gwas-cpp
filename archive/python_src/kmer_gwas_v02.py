#!/usr/bin/env python
#
#
#from multiprocessing import Pool
from multiprocessing.pool import ThreadPool as Pool
import math
import pprint, psutil, os
from p_tqdm import p_map
import time


pp = pprint.PrettyPrinter(indent=4)
process = psutil.Process(os.getpid())

k=71
CGR_SIZE = 32
NPROC = 10 # expose to user
kmers = dict()

myFile="test_1m.fq"
#myFile="test_1m.fq"
#myFile="/local/data/kaust_research/adil_salhi/kmer_gwas/brande_wulff_gwas/AEG0393_1.fq"
## CGR dicts
X={'A':0, 'C':0, 'G': 1, 'T': 1}
Y={'A':0, 'C':1, 'G': 1, 'T': 0}
rc_dict = {'A':'T','T':'A','G':'C','C':'G'}

def canonical(seq):
	rc = ''.join(list(map(rc_dict.get, seq[::-1])))
	return seq if seq < rc else rc

##############################################
def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n
##############################################
def get_kmers (seq):

	for i in range(len(seq) - k + 1):
		kmer = seq[i:i+k]
		if "N" in kmer:
			continue
		kmer = canonical(kmer)
		kmers[kmer_to_cgr(kmer)] = kmers.get(kmer, 0) + 1
		#kmers[kmer] = kmers.get(kmer, 0)+1
##############################################

##############################################
## calculate CGR coordinates 
def kmer_to_cgr (kmer):
	x, y = (0, 0)
	for i in range(len(kmer)):
		if i == 0:
			x = X[kmer[i]]
			y = Y[kmer[i]]
		else:
			x = (x + X[kmer[i]])/2
			y = (y + Y[kmer[i]])/2
	x = int(x * math.pow(2, CGR_SIZE))
	y = int(y * math.pow(2, CGR_SIZE))

	key = int(x + y * math.pow(2, CGR_SIZE))
	return  key
##############################################
##
#def process_file_parallel ():


##############################################
print("##############################################")
tic = time.perf_counter()
cts = time.process_time()
print("reading "+myFile+" ...", flush=True)
with open (myFile, 'r+') as f:
	lines= f.readlines()
toc = time.perf_counter()
cte = time.process_time()
print(f"Read {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
print(f"Read {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)
##############################################
## AEG0393_1.fq  --> 58gb
## RSS - After Reading File: 103.0G
##############################################
print("RSS - After Reading File: "+str(bytes2human(process.memory_info().rss)), flush=True)
#print('RAM after reading file % used:', psutil.virtual_memory()[2], flush=True)
pp.pprint(psutil.virtual_memory())
print("##############################################")
# Get the sequence lines only from fastq
print("keep sequences data only ...")
##############################################
## AEG0393_1.fq  --> 58gb
## RSS - After Pruning File: 35.0G 
## 100.41user 135.89system 3:56.32elapsed 99%CPU (0avgtext+0avgdata 145947316maxresident)k
tic = time.perf_counter()
cts = time.process_time()
lines=list(map(str.strip, lines[1::4]))
toc = time.perf_counter()
cte = time.process_time()
print(f"Pruned {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
print(f"Pruned {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)
##############################################
## AEG0393_1.fq  --> 58gb
## RSS - After Pruning File: 35.0G 
## 152.13user 153.68system 5:04.89elapsed 100%CPU (0avgtext+0avgdata 146735768maxresident)k 
#with Pool() as p:
#	lines=list(p.map(str.strip, lines[1::4]))
##############################################
## AEG0393_1.fq  --> 58gb
## RSS - After Pruning File: 35.0G 
#lines=p_map(str.strip, lines[1::4])
##############################################
print("RSS - After Pruning File: "+str(bytes2human(process.memory_info().rss)), flush=True)
#print('RAM after % used:', psutil.virtual_memory()[2], flush=True)
pp.pprint(psutil.virtual_memory())
print("##############################################")
# _=list(map(get_kmers, lines))
print("get kmers ...")
# with Pool(NPROC) as p:
print("#######################")
### p_map
print("### p_map")
tic = time.perf_counter()
cts = time.process_time()
_=p_map(get_kmers, lines, num_cpus=NPROC)
#with Pool(NPROC) as p:
#	_=list(p.map(get_kmers, lines))
toc = time.perf_counter()
cte = time.process_time()
print(f"Calculated Kmers for {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
print(f"Calculated Kmers for {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)

print("RSS - After Getting Kmers: "+str(bytes2human(process.memory_info().rss)), flush=True)
#print('RAM after kmers % used:', psutil.virtual_memory()[2], flush=True)
pp.pprint(psutil.virtual_memory())
print("#######################")
### Pool
print("### Pool")
kmer = {}
tic = time.perf_counter()
cts = time.process_time()
#_=p_map(get_kmers, lines, num_cpus=NPROC)
with Pool(NPROC) as p:
	_=list(p.map(get_kmers, lines))
toc = time.perf_counter()
cte = time.process_time()
print(f"Calculated Kmers for {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
print(f"Calculated Kmers for {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)
print("RSS - After Getting Kmers: "+str(bytes2human(process.memory_info().rss)), flush=True)
#print('RAM after kmers % used:', psutil.virtual_memory()[2], flush=True)
pp.pprint(psutil.virtual_memory())
print("#######################")
#map(get_kmers, lines)
print("printing kmers ...DONE", flush=True)
#pp.pprint(kmers)
print("##############################################")



	
