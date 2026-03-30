#!/usr/bin/env python
#
#
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool as Threadpool
import math
import pprint, psutil, os, sys
from p_tqdm import p_map
from tqdm import tqdm
import time
import gc
#from collections import Counter

gc.disable()

pp = pprint.PrettyPrinter(indent=4)
process = psutil.Process(os.getpid())

k=71
CGR_SIZE = 32
NPROC = 30 # expose to user
NUM_CHUNKS=120000
kmers = dict()

#myFile="test_2m.fq"
#myFiles=["test_3m.fq","test_4m.fq","test_5m.fq"]
#myFiles=["test_3m.fq"]
myFiles=["/local/data/kaust_research/adil_salhi/kmer_gwas/brande_wulff_gwas/AEG0393_1.fq"]
#myFiles=["test_10m.fq"]
# myFiles=["test_1m.fq"]
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
def get_kmers (chunk_record):
	local_kmer_bag = dict()
	global kmer_bag
	rec, chunk = chunk_record
	# for seq in chunk:
	for a in range(len(chunk)):
		for i in range(len(chunk[a]) - k + 1):
			kmer = chunk[a][i:i+k]
			if "N" in kmer:
				continue
			kmer = canonical(kmer)
			#local_kmer_bag[kmer_to_cgr(kmer)] = kmers.get(kmer, 0) + 1
			local_kmer_bag[kmer] = kmers.get(kmer, 0) + 1
		#kmers[kmer] = kmers.get(kmer, 0)+1
		chunk[a] = 0 # release memory
	# print(f"DONE chunk #{rec}", flush=True)
	return local_kmer_bag
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

##############################################
def update_global_dict (current_dict):
	global global_kmer_index
	for key in current_dict:
		global_kmer_index[key] = global_kmer_index.get(key, 0) + current_dict[key]

##############################################


##############################################
for myFile in myFiles:
	print(f"*********** START {myFile} **************************")
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
	print("### multiprocessing")
	tic = time.perf_counter()
	cts = time.process_time()

	chunk_size = int(len(lines)/NUM_CHUNKS)+1
	print(f"chunk size: {chunk_size} NPROC {NPROC}")
	kmer_bag = list()
	data=list()
	for i in range(NUM_CHUNKS):
		data.append(lines[i*chunk_size:(i+1)*chunk_size])
	del lines
	data=list(enumerate(data))
	# kmer_bag = p_map(get_kmers, data,  num_cpus=NPROC)
	with Pool(NPROC) as p:
		kmer_bag = list(p.map(get_kmers, data))

	toc = time.perf_counter()
	cte = time.process_time()
	print(f"Calculated Kmers for {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
	print(f"Calculated Kmers for {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)

	print("RSS - After Getting Kmers: "+str(bytes2human(process.memory_info().rss)), flush=True)
	#print('RAM after kmers % used:', psutil.virtual_memory()[2], flush=True)
	pp.pprint(psutil.virtual_memory())
	print("#######################")

	print("##############################################")
	print("Merging chunks ...")
	tic = time.perf_counter()
	cts = time.process_time()

	### REFERENCE ###################
	global_kmer_index=dict()
	for i in tqdm(range(len(kmer_bag))):
		# print(f"adding dict {i} with {len(kmer_bag[i].keys())} keys")
		for key in kmer_bag[i]:
			global_kmer_index[key] = global_kmer_index.get(key, 0) + kmer_bag[i][key]
		kmer_bag[i] = dict()
	###############################

	toc = time.perf_counter()
	cte = time.process_time()
	print(f"Calculated Kmers for {myFile} in {toc - tic:0.4f} seconds (Wallclock)", flush=True)
	print(f"Calculated Kmers for {myFile} in {cts - cte:0.4f} seconds (CPU Time)", flush=True)
	print("RSS - After Getting Kmers: "+str(bytes2human(process.memory_info().rss)), flush=True)


	print(f"number of kmers {len(global_kmer_index.keys())}")
	#pp.pprint(global_kmer_index)
	print("##############################################")

	print(f"*********** END {myFile} **************************")
