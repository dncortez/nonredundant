import os
import time
import sys
import numpy as np

start_time=time.time()

arguments=sys.argv
options={arguments[i]:(arguments[i+1] if i<len(arguments)-1 and arguments[i+1][:2]!="--" and arguments[i] in ["--fmt"] else None) for i in range(len(arguments)) if arguments[i][:2]=="--"}
non_kw=[arg for arg in arguments[1:] if not (arg in options.values() or arg[:2]=="--")]

if "--help" in options:
	script_dir=os.path.dirname(__file__)+"/"
	help_doc=open(script_dir+"help.txt").readlines()
	exit("".join(help_doc))

ext=("."+options["--fmt"]) if "--fmt" in options and not options["--fmt"] is None else ".faa"

folder=os.getcwd()+"/"+(non_kw[0] if len(non_kw)>0 else "")+"/"

if not os.path.exists(folder): exit("Input folder "+str(folder[:-1])+" does not exist")

out_folder=folder+"nr/"
out_nrfile=out_folder+"nr.faa"

def readfasta(filepath, mode=0, fullid=False):
	inp=open(filepath,"r").readlines()
	out, outids = [], []
	for line in inp:
		if line[0]==">":
			if mode in (0, 2): out.append("")
			if mode in (1, 2): outids.append(line[1:-1].split(" ")[0])
		elif  mode in (0, 2):
			out[-1]+=line[:-1]
	return out if mode==0 else (outids if mode==1 else (out, outids))
def splitbylen(strings):
	lengths=list(set(sorted([len(string) for string in strings])))
	out, idxout= [[] for i in lengths], [[] for i in lengths]
	for i, string in enumerate(strings):
		len_ind=lengths.index(len(string))
		out[len_ind].append(string)
		idxout[len_ind].append(i)
	return out, idxout
def donr(folder, enforce=False, critical_mass=100, check=True, verbose=True, silent=False):
	if enforce or not os.path.exists(out_nrfile):
		print("Reading proteomes...")
		filenames=[file[:-len(ext)] for file in os.listdir(folder) if file[-len(ext):]==ext]
		if len(filenames)==0: exit("No proteomes in target folder. Proteome files must have a .faa extension. Try specifying a path or use the --fmt option. See --help for more details")
		if not os.path.exists(out_folder): os.makedirs(out_folder)
		ids, proteomes= [], []
		for file_i, filename in enumerate(filenames):
			readprot=readfasta(folder+filename+ext)
			proteomes=proteomes+readprot
			ids+=[(file_i, prot_i) for prot_i in range(len(readprot))]
			print(file_i, "proteomes read so far") if file_i%( max([5]+[f for f in [10,20,50,100] if len(filenames)/f>4]) )==0 and verbose else None
		
		
		nrids=[-1]*len(ids)
		
		if verbose: print("Number of read proteomes:", len(filenames), ". Total proteins:",len(proteomes),"\nSpliting proteomes by protein length")
		
		proteomeses, indexes = splitbylen(proteomes)
		sizes=[len(proteome) for proteome in proteomeses]
		
		if verbose: print("Proteomes split into "+str(len(proteomeses))+" clusters of same length proteins. Biggest cluster: "+str(max(sizes))+" proteins of length "+str(len(proteomeses[sizes.index(max(sizes))][0]))+" aminoacids")

		nrout=open(out_nrfile,"w")
		
		def nrsamelength(proteins, index, snrid):
			nr, index, proteins = [], [index], [proteins]
			listo=False
			l, lmax = 0, len(proteins[0][0])
			while l<lmax and not listo:
				nproteins, nindex = [], []
				for prots, idxs in zip(proteins, index):
					if len(prots)==1:
						nr.append(prots[0])
						nrids[idxs[0]]=snrid
						snrid+=1
					elif len(prots)<=critical_mass:
						nrs, nrsidxs = [], []
						for prot, id in zip(prots, idxs):
							if not prot in nrs:
								nrs.append(prot)
								nrsidxs.append(snrid)
								nrids[id]=snrid
								snrid+=1
							else:
								nrids[id]=nrsidxs[nrs.index(prot)]
						nr=nr+nrs
					else:
						nprots, nidxs, letras = [], [], ""
						for prot, id in zip(prots, idxs):
							if prot[l] in letras:
								nprots[letras.index(prot[l])].append(prot)
								nidxs[letras.index(prot[l])].append(id)
							else:
								letras+=prot[l]
								nprots.append([prot])
								nidxs.append([id])
						nproteins=nproteins+nprots
						nindex=nindex+nidxs
				proteins=nproteins
				index=nindex
				if len(proteins)==0:
					listo=True
				l+=1
			if len(proteins)>0:
				for p in range(0,len(proteins)):
					prots=proteins[p]
					nr.append(prots[0])
					for i in index[p]:
						nrids[i]=snrid
					snrid+=1
			return nr
		
		startnrid, prnrid, porc, acum = 0,0,0,0
		for prots, size, indexs in zip(proteomeses, sizes, indexes):
			acum+=size
			if verbose and not int(acum*20.0/sum(sizes))==porc:
				porc=int(acum*20.0/sum(sizes))
				print("Creating Non Redundant proteome... "+str(porc*5)+"%")
			nr=nrsamelength(prots, indexs, startnrid)
			startnrid+=len(nr)

			for p in range(0,len(nr)):
				nrout.write(">NR_000000000"[:-len(str(prnrid))]+str(prnrid)+"\n")
				prnrid+=1
				for a in range(0,len(nr[p]),80):
					nrout.write(nr[p][a:a+80]+"\n")
		nrout.close()
		
		print("Non redundant proteome made. Total proteins in original proteomes (redundant): "+str(len(proteomes))+". Non redundant proteome size: "+str(prnrid))
		if verbose: print("Creating total proteins table")

		protable=[0]*len(nrids)
		for i in range(len(nrids)): protable[i]=[filenames[ids[i][0]],ids[i][1],"",nrids[i]]
		
		if verbose: print("Pairing proteins with NR id")
		nprot=0
		for filename in filenames:
			protids=readfasta(folder+filename+ext, mode=1)
			for l in range(len(protids)): protable[nprot+l][2]=protids[l]
			nprot+=len(protids)
		
		if verbose: print("Exporting table")
		
		nridsout=open(out_folder+"nr_table.tsv","w")
		nridsout.write("proteome\tn_protein\tprotein_id\tnr_id\n")
		for line in protable:
			nridsout.write(str(line[0])+"\t"+ str(line[1])+"\t"+line[2]+"\t"+str(line[3])+"\n")
		nridsout.close()
		
		if check:
			checkredundancy(verbose=verbose)
			checkcorrespondence(verbose=verbose)
	else: print(out_folder+" already exists. Try changing the target directory or using the --enforce option")
def checkredundancy(verbose=True):
	print("Checking redundancy")
	nr=readfasta(out_nrfile)
	sizebookmark=0
	psize=0
	for i in range(len(nr)):
		if not len(nr[i])==psize:
			psize=len(nr[i])
			sizebookmark=i
			if verbose: print("Analyzing proteins of size "+str(psize))
		prange=nr[sizebookmark:i]
		if nr[i] in prange:
			print("Redundancy found! Protein "+str(i)+" is identical to protein "+str(prange.index(nr[i])))
			exit()
	print("The NR Database is non redundant!")
def checkcorrespondence(verbose=True):
	print("Checking that every protein in the organisms' proteomes are in nr file.")
	if verbose: print("Loading nr table")
	nr_table=[line[:-1].split("\t") for line in open(out_folder+"nr_table.tsv").readlines()[1:]]
	nr=readfasta(out_nrfile)
	pfol=""
	proteome=[]
	pindex=0
	protsanalyzed=0
	totalprots=0
	for line in nr_table:
		if not line[0]==pfol:
			if len(proteome)>0:
				if pindex+1==len(proteome):
					if verbose: print("Proteome "+str(pfol)+" OK")
				else:
					print("Proteome "+str(pfol)+" has "+str(len(proteome))+" but "+str(pindex+1)+" proteins were analyzed")
			pfol=line[0]
			proteome=readfasta(folder+line[0]+ext)
			totalprots+=len(proteome)
		pindex=int(line[1])
		nrindex=int(line[3])
		if not proteome[pindex]==nr[nrindex]:
			print("ERROR!!! Protein "+str(pindex)+" of genome "+str(pfol)+" (ID: "+str(line[2])+") differs to protein "+str(nrindex)+" of NR:")
			print(">"+str(line[2])+" "+str(pfol)+"_p"+"000000000"[:-len(str(pindex))]+str(pindex))
			print(proteome[pindex])
			print(">NR_"+"000000000"[:-len(str(nrindex))]+str(nrindex))
			print(nr[nrindex])
			exit()
		protsanalyzed+=1
	print("Total proteins in genomes/Proteins with a NR analog: "+str(totalprots)+"/"+str(protsanalyzed)+". All proteins are in the NR database")
	
if "--justckeck" in options:
	checkredundancy(verbose=not "--silent" in options)
	checkcorrespondence(verbose=not "--silent" in options)
else:
	donr(folder, enforce=True, critical_mass=100, check="--check" in options, verbose=not "--silent" in options)


seconds=time.time()-start_time
if not "--silent" in options: print("Done. Elapsed time:", seconds//3600, "hours,", (seconds//60)%60, "minutes,", int(seconds)%60, "seconds and", int((seconds%1)*1000), "milliseconds.\n:)")