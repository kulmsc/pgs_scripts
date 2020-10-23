import numpy as np
import datetime
import time
import pickle
import pdb
import gzip
import os
import sys
import sys
import re

author = "Bentham"

def normRead(fileName, withHeader = True, delim = '\t', removeQuote = False):
	with open(fileName,"r",encoding = "latin-1") as f:
		totalData=[]
		for line in f.read().splitlines():
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

date_ass,ass_header = normRead("date_assessed.csv", True, ",", True)
meds, meds_header = normRead("medications.csv", True, ",", True)
self_rep, self_rep_header = normRead("self_report_diag.csv", True, ",", True)
cancer_rep, cancer_rep_header = normRead("self_report_cancer.csv", True, ",", True)
big_eid, eid_head = normRead("eid.csv", True, ",", True)

date_ass = date_ass[big_eid[:,0].argsort(),:]
meds = meds[big_eid[:,0].argsort(),:]
self_rep = self_rep[big_eid[:,0].argsort(),:]
cancer_rep = cancer_rep[big_eid[:,0].argsort(),:]
big_eid = big_eid[big_eid[:,0].argsort(),:]
defs, head_defs = normRead("../descript_defs/author_defs")

cancer_phase = [x.split("-")[1].split(".")[0] for x in cancer_rep_header]
self_phase = [x.split("-")[1].split(".")[0] for x in self_rep_header]
meds_phase = [x.split("-")[1].split(".")[0] for x in meds_header]

diag, head_diag = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt")
oper, head_oper = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_oper.txt")
hesin, head_hesin = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin.txt")
diag = diag[diag[:,0].argsort(),:]
oper = oper[oper[:,0].argsort(),:]
hesin = hesin[hesin[:,0].argsort(),:]
diag_eid = np.unique(diag[:,0])
oper_eid = np.unique(diag[:,0])
diag_max = diag.shape[0] - 1
oper_max = oper.shape[0] - 1
diag_eid_list = diag[:,0].tolist()
oper_eid_list = oper[:,0].tolist()

df_occur = np.zeros((big_eid.shape[0], 6))
df_date = np.tile("_________", (big_eid.shape[0], 6))


def_ind = defs[:,0].tolist().index(author)
split_defs = list()
for type_def in defs[def_ind, 3:9]:
	split_defs.append(type_def.split("|"))
use_defs = dict(zip(head_defs[3:9], split_defs))

cancer_sr_time = 0
noncancer_sr_time = 0
hesin_setup_time = 0
icd9_time = 0
icd10_time = 0
hesin_oper_time = 0
med_time = 0

start_eid_diag = 0
start_eid_oper = 0

for ind in range(big_eid.shape[0]):
	if ind % 10000 == 0:
		print(ind)

	curr_eid = big_eid[ind][0]

	#CANCER SELF REPORT
	if use_defs["cancer"][0] != "NA":
		t1 = time.time()
		if any(np.isin(use_defs["cancer"], cancer_rep[ind,:])):
			locs = np.where(np.isin(cancer_rep[ind,:], use_defs["cancer"]))[0]
			df_date[ind,0] = date_ass[ind, int(cancer_phase[locs[0]])]
			df_occur[ind,0] = 1

		t2 = time.time()
		cancer_sr_time += t2-t1

	#NON-CANCER SELF REPORT
	if use_defs["noncancer"][0] != "NA":
		t22 = time.time()
		if any(np.isin(use_defs["noncancer"], self_rep[ind,:])):
			locs = np.where(np.isin(self_rep[ind,:], use_defs["noncancer"]))[0]
			df_date[ind,1] = date_ass[ind, int(self_phase[locs[0]])]
			df_occur[ind,1] = 1

		t3 = time.time()
		noncancer_sr_time += t3-t22


	#HESIN DIAG
	if curr_eid in diag_eid_list:
		t4 = time.time()
		sub_hesin = hesin[hesin[:,0] == curr_eid,:]
		if curr_eid != diag_eid_list[-1]:
			next_eid = np.unique(diag[start_eid_diag:(start_eid_diag+7000),0])[1]
			ext_eid = diag_eid_list.index(next_eid)
		else:
			ext_eid = diag_max

		t5 = time.time()
		hesin_setup_time += t5-t4

		#ICD - 9
		use_short_icd9_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,4]]
		if any(np.isin(use_defs["icd9"], use_short_icd9_diag)):
			short_icd9_locs = np.where(np.isin(use_short_icd9_diag, use_defs["icd9"]))[0][0]
		else:
			short_icd9_locs = 10000
		if any(np.isin(use_defs["icd9"], diag[start_eid_diag:ext_eid,4])):
                        long_icd9_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,4], use_defs["icd9"]))[0][0]
		else:
			long_icd9_locs = 10000
		if long_icd9_locs != 10000 or short_icd9_locs != 10000:
			icd9_locs = min((long_icd9_locs, short_icd9_locs))
			icd9_ins_index = diag[start_eid_diag+icd9_locs,1]

			raw_date = sub_hesin[sub_hesin[:,1] == icd9_ins_index, 5][0]
			if raw_date == "":
				raw_date = sub_hesin[sub_hesin[:,1] == icd9_ins_index, 21][0]
				if raw_date == "":
					raw_date = sub_hesin[sub_hesin[:,1] == icd9_ins_index, 4][0]

			df_occur[ind,2] = 1
			df_date[ind,2] = raw_date
		
		t6 = time.time()
		icd9_time += t6-t5

		#ICD - 10
		use_short_icd10_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,6]]
		if any(np.isin(use_defs["icd10"], use_short_icd10_diag)):
			short_icd10_locs = np.where(np.isin(use_short_icd10_diag, use_defs["icd10"]))[0][0]
		else:
			short_icd10_locs = 10000
		if any(np.isin(use_defs["icd10"], diag[start_eid_diag:ext_eid,6])):
			long_icd10_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,6], use_defs["icd10"]))[0][0]
		else:
			long_icd10_locs = 10000
		if long_icd10_locs != 10000 or short_icd10_locs != 10000:
			icd10_locs = min((long_icd10_locs, short_icd10_locs))
			icd10_ins_index = diag[start_eid_diag+icd10_locs,1]

			raw_date = sub_hesin[sub_hesin[:,1] == icd10_ins_index, 5][0]
			if raw_date == "":
				raw_date = sub_hesin[sub_hesin[:,1] == icd10_ins_index, 21][0]
				if raw_date == "":
					raw_date = sub_hesin[sub_hesin[:,1] == icd10_ins_index, 4][0]

			df_occur[ind,3] = 1
			df_date[ind,3] = raw_date
		start_eid_diag = ext_eid

		t7 = time.time()
		icd10_time += t7-t6


	#HESIN OPER
	if curr_eid in oper_eid_list and use_defs["opcs"][0] != "NA":
		t8 = time.time()
		if curr_eid != oper_eid_list[-1]:
			next_eid = np.unique(diag[start_eid_oper:(start_eid_oper+7000),0])[1]
			ext_eid = oper_eid_list.index(next_eid)
		else:
			ext_eid = oper_max

		use_short_oper = [x[0:3] for x in oper[start_eid_oper:ext_eid,7]]
		if any(np.isin(use_defs["opcs"], use_short_oper)):
			short_oper_locs = np.where(np.isin(use_short_oper, use_defs["opcs"]))[0][0]
		else:
			short_oper_locs = 10000
		if any(np.isin(use_defs["opcs"], oper[start_eid_oper:ext_eid,7])):
			long_oper_locs = np.where(np.isin(oper[start_eid_oper:ext_eid,7], use_defs["opcs"]))[0][0]
		else:
			long_oper_locs = 10000
		if long_oper_locs != 10000 or short_oper_locs != 10000:
			oper_locs = min((long_oper_locs, short_oper_locs))
			oper_ins_index = oper[start_eid_oper+oper_locs,1]

			raw_date = sub_hesin[sub_hesin[:,1] == oper_ins_index, 5][0]
			if raw_date == "":
				raw_date = sub_hesin[sub_hesin[:,1] == oper_ins_index, 21][0]
				if raw_date == "":
					raw_date = sub_hesin[sub_hesin[:,1] == oper_ins_index, 4][0]

			df_occur[ind,4] = 1
			df_date[ind,4] = raw_date
		start_eid_oper = ext_eid

		t9 = time.time()
		hesin_oper_time += t9-t8

	#MEDICATION
	if use_defs["meds"][0] != "NA":
		t10 = time.time()
		if any(np.isin(use_defs["meds"], meds[ind,:])):
			locs = np.where(np.isin(meds[ind,:], use_defs["opcs"]))[0]
			df_date[ind,5] = date_ass[ind, int(meds_phase[locs[0]])]
			df_occur[ind,5] = 1
		t11 = time.time()
		med_time += t11 - t10



np.savetxt("diag.coding.txt.gz", df_occur, fmt='%i')
np.savetxt("diag.time.txt.gz", df_date, fmt='%s')



print("end")
