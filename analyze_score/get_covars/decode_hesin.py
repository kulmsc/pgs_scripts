import numpy as np
import datetime
import time
import pickle
import pdb
import gzip
import os
import sys
import re
import gzip

author = sys.argv[1]
#author = "Christophersen"
phen_method = "all"          
#Still need to sort everything so that the extEid process works

def normRead(fileName, withHeader = True, delim = '\t', removeQuote = False):
	totalData = []
	if fileName[-2:] == "gz":
		with gzip.open(fileName,"r") as f:
			for line in f:
				totalData.append(line.decode().strip().split(delim))

	else:
		with open(fileName,"r",encoding = "latin-1") as f:
			for line in f.read().splitlines():
				#for line in f:
				if removeQuote:
					line = line.replace('"', '').strip()
				totalData.append(line.split(delim))
	if withHeader:
		header=totalData[0]
		del totalData[0]
	else:
		header = None
	totalData=np.array(totalData)
	return(totalData,header)


author_needs, no_header = normRead("../descript_defs/author_to_covar_hesin", False)
author_needs = author_needs[author_needs[:,0] == author, 1]
author_needs = author_needs[0].split(",")
if author_needs[0] == "NA":
	print("no covars")
	sys.exit()

#pdb.set_trace()
covar_defs, covar_header = normRead("../descript_defs/covar_defs_hesin", False)
use_defs = []
for i in range(len(author_needs)):
	split_defs = [j.split(",") for j in covar_defs[covar_defs[:,0] == author_needs[i],2:][0]]
	use_defs.append(dict(zip(["ICD10", "ICD9", "selfrep", "cancer"], split_defs)))

#for i in range(covar_defs.shape[0]):
#	split_defs = [j.split(",") for j in covar_defs[i,2:]]
#	use_defs.append(dict(zip(["ICD10", "ICD9", "selfrep", "cancer"], split_defs)))


date_ass,ass_header = normRead("../construct_defs/date_assessed.csv", True, ",", True)
self_rep, self_rep_header = normRead("../construct_defs/self_report_diag.csv", True, ",", True)
cancer_rep, cancer_rep_header = normRead("../construct_defs/self_report_cancer.csv", True, ",", True)
big_eid, eid_head = normRead("../construct_defs/eid.csv", True, ",", True)

date_ass = date_ass[big_eid[:,0].argsort(),:]
self_rep = self_rep[big_eid[:,0].argsort(),:]
cancer_rep = cancer_rep[big_eid[:,0].argsort(),:]
big_eid = big_eid[big_eid[:,0].argsort(),:]

cancer_phase = [x.split("-")[1].split(".")[0] for x in cancer_rep_header]
self_phase = [x.split("-")[1].split(".")[0] for x in self_rep_header]

diag, head_diag = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt")
hesin, head_hesin = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin.txt")
diag = diag[diag[:,0].argsort(),:]
hesin = hesin[hesin[:,0].argsort(),:]
diag_eid = np.unique(diag[:,0])
diag_max = diag.shape[0] - 1
diag_eid_list = diag[:,0].tolist()
uni_diag_eid_list = [x for x in set(diag_eid_list)]


print("reading")
known_diag, known_head = normRead("../construct_defs/pheno_defs/diag." + author.lower() + ".txt.gz", False, " ")
known_time, known_head = normRead("../construct_defs/pheno_defs/time." + author.lower() + ".txt.gz", False, " ")
new_time = np.tile(datetime.datetime.strptime("31/12/2020", '%d/%m/%Y'), (known_time.shape[0], known_time.shape[1]))
print("read in")
for i in range(known_time.shape[0]):
	for j in range(known_time.shape[1]):
		if known_time[i,j] != '__________':
			if j == 0 or j == 1 or j == 5:
				new_time[i,j] = datetime.datetime.strptime(known_time[i,j], "%Y-%m-%d")
			else:
				new_time[i,j] = datetime.datetime.strptime(known_time[i,j], "%d/%m/%Y")
known_time = new_time
	
if phen_method == "icd":
	known_diag = know_diag[:,2:4]
	known_time = known_time[:,2:4]
elif phen_method == "selfrep":
	known_diag = known_diag[:,0:2]
	known_time = known_time[:,0:2]
elif phen_method == "icd_selfrep":
	known_diag = known_diag[:,0:4]
	known_time = known_time[:,0:4]
elif phen_method == "double":
	known_diag = known_diag.astype("int")
	for_double = np.sum(known_diag, axis = 1)

final_diag = []
final_time = []
for i in range(known_diag.shape[0]):
	if phen_method == "double":
		if for_double[i] > 1:
			final_diag.append(1)
			final_time.append(min(known_time[i,:]))
		else:
			final_diag.append(0)
			final_time.append("NA")
	else:
		if "1" in known_diag[i,:]:
			final_diag.append(1)
			final_time.append(min(known_time[i,:]))
		else:
			final_diag.append(0)
			final_time.append("NA")


#STILL NEED TO FIGURE OUT HOW I WANT TO SPLIT THIS - right now assume I dont

start_eid_diag = 0
df_occur = np.zeros((big_eid.shape[0], len(use_defs)))
df_date = np.tile("_________", (big_eid.shape[0], len(use_defs)))

hesin_0_col = hesin[:,0].astype("int")
hesin_1_col = hesin[:,1].astype("int")
diag_1_col = diag[:,1].astype("int")


cancer_time = 0
selfrep_time = 0
setuphesin_time = 0
icd10_time = 0
icd9_time = 0
icd10_1_time = 0
icd10_2_time = 0
icd10_3_time = 0

for ind in range(big_eid.shape[0]):
	#print(ind)
	get_ext_eid = True
	curr_eid = big_eid[ind][0]
	if ind % 1000 == 0:
		print(ind)
		#pdb.set_trace()

	for trait_ind in range(len(use_defs)):
		#print("trait_ind", trait_ind)

		#t1 = time.time()
		#CANCER SELF REPORT
		if use_defs[trait_ind]["cancer"] != "NA":
			if any(np.isin(use_defs[trait_ind]["cancer"], cancer_rep[ind,:])):
				locs = np.where(np.isin(cancer_rep[ind,:], use_defs[trait_ind]["cancer"]))[0]
				poss_date = date_ass[ind, int(cancer_phase[locs[0]])]
				if poss_date == "":
					poss_date = date_ass[ind, 0]
				poss_date = datetime.datetime.strptime(poss_date, "%Y-%m-%d")
				if final_time[i] == "NA":
					df_occur[ind, trait_ind] = 1
					df_date[ind, trait_ind] = poss_date
				elif poss_date < final_time[i]:
					df_occur[ind, trait_ind] = 1
					df_date[ind, trait_ind] = poss_date

		#t2 = time.time()
		#cancer_time += t2-t1

		#NON-CANCER SELF REPORT
		if use_defs[trait_ind]["selfrep"] != "NA":
			if any(np.isin(use_defs[trait_ind]["selfrep"], self_rep[ind,:])):
				locs = np.where(np.isin(self_rep[ind,:], use_defs[trait_ind]["selfrep"]))[0]
				poss_date = date_ass[ind, int(self_phase[locs[0]])]

				if poss_date == "":
					poss_date = date_ass[ind, 0] #problem at in 7009
				poss_date = datetime.datetime.strptime(poss_date, "%Y-%m-%d")
				if final_time[i] == "NA":
					df_occur[ind, trait_ind] = 1
					df_date[ind, trait_ind] = poss_date
				elif poss_date < final_time[i]:
					df_occur[ind, trait_ind] = 1
					df_date[ind, trait_ind] = poss_date

		#t3 = time.time()
		#selfrep_time += t3-t2


		#HESIN DIAG
		if curr_eid in uni_diag_eid_list:
			#pdb.set_trace()
			#t4 = time.time()

			if get_ext_eid:
				if curr_eid != diag_eid_list[-1]:
					next_eid = np.unique(diag[start_eid_diag:(start_eid_diag+7000),0])[1]
					ext_eid = diag_eid_list.index(next_eid)
				else:
					ext_eid = diag_max
				get_ext_eid = False

			#t5 = time.time()
			#setuphesin_time += t5-t4

			#ICD -9
			use_short_icd9_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,4]]
			if any(np.isin(use_defs[trait_ind]["ICD9"], use_short_icd9_diag)):
				short_icd9_locs = np.where(np.isin(use_short_icd9_diag, use_defs[trait_ind]["ICD9"]))[0][0]
			else:
				short_icd9_locs = 10000
			if any(np.isin(use_defs[trait_ind]["ICD9"], diag[start_eid_diag:ext_eid,4])):
				long_icd9_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,4], use_defs[trait_ind]["ICD9"]))[0][0]
			else:
				long_icd9_locs = 10000
			if long_icd9_locs != 10000 or short_icd9_locs != 10000:
				icd9_locs = min((long_icd9_locs, short_icd9_locs))
				icd9_ins_index = diag[start_eid_diag+icd9_locs,1]

				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),5][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),21][0]
					if raw_date == "":
						raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),4][0]

				df_occur[ind, trait_ind] = 1
				df_date[ind, trait_ind] = raw_date

			#t6 = time.time()
			#icd9_time += t6 - t5


			#ICD - 10
			use_short_icd10_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,6]]
			if any(np.isin(use_defs[trait_ind]["ICD10"], use_short_icd10_diag)):
				short_icd10_locs = np.where(np.isin(use_short_icd10_diag, use_defs[trait_ind]["ICD10"]))[0][0]
			else:
				short_icd10_locs = 10000
			#t61 = time.time()
			#icd10_1_time += t61 - t6
			if any(np.isin(use_defs[trait_ind]["ICD10"], diag[start_eid_diag:ext_eid,6])):
				long_icd10_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,6], use_defs[trait_ind]["ICD10"]))[0][0]
			else:
				long_icd10_locs = 10000
			#t62 =time.time()
			#icd10_2_time += t62 - t61
			if long_icd10_locs != 10000 or short_icd10_locs != 10000:
				#pdb.set_trace()
				t63 = time.time()
				icd10_locs = min((long_icd10_locs, short_icd10_locs))
				#icd10_ins_index = diag[start_eid_diag+icd10_locs,1]
				#icd10_ins_index = diag_1_col[start_eid_diag+icd10_locs]

				#raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),5][0]
				raw_date = hesin[np.logical_and(hesin_0_col == int(curr_eid), hesin_1_col == diag_1_col[start_eid_diag+icd10_locs]),5][0]
				#raw_date = hesin[~((hesin_0_col == int(curr_eid)) & (hesin_1_col == diag_1_col[start_eid_diag+icd10_locs])),5][0]
				#t64 = time.time()
				#icd10_3_time += t64 - t63
				if raw_date == "":
					#raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin_1_col == diag_1_col[start_eid_diag+icd10_locs]),21][0]
					raw_date = hesin[np.logical_and(hesin_0_col == int(curr_eid), hesin_1_col == diag_1_col[start_eid_diag+icd10_locs]),21][0]
					#raw_date = hesin[~((hesin_0_col == int(curr_eid)) & (hesin_1_col == diag_1_col[start_eid_diag+icd10_locs])),21][0]
					if raw_date == "":
						#raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin_1_col == diag_1_col[start_eid_diag+icd10_locs]),4][0]
						raw_date = hesin[np.logical_and(hesin_0_col == int(curr_eid), hesin_1_col == icd10_ins_index),4][0]
						#raw_date = hesin[~((hesin_0_col == int(curr_eid)) & (hesin_1_col == diag_1_col[start_eid_diag+icd10_locs])),4][0]

				df_occur[ind, trait_ind] = 1
				df_date[ind, trait_ind] = raw_date

			#t7 = time.time()
			#icd10_time += t7 - t6


		start_eid_diag = ext_eid

#save the output phenos
np.savetxt("hesin_covars/diag.coding." + author.lower() + ".txt.gz", df_occur, fmt='%i')
np.savetxt("hesin_covars/diag.time." + author.lower() + ".txt.gz", df_date, fmt='%s')








