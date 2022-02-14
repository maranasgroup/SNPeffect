import numpy as np 
import os 
import sys
from os import listdir 

sys.path.append('.')

###define useful paths 
data_dir = './../data/'
###folder where the metabolic model is stored
gsm_dir = data_dir+'gsm/'
##directory where genotype-wise metabolite levels are stored
metabdatapath = data_dir+'metabolomics/' 
###directory where genotype-wise sequence data is stored
snpdatapath = data_dir+'genomics/'
##directory where the output files are stored
respath = './../out/'
respath_forgams = 'gams_inputs/'


'''
this code needs the following inputs 


1) data_dir
>> AllGenotypesInStudy.txt - list of genotypes in the study
>> Biomass_Cons_AllGenotypes.txt - relative growth rate of each genotype wrt the reference genotype
>> Biomass_ConsList_AllGenotypes.txt - list of constraints for the above file


2) gsm_dir
>> smbl_metabolites.txt - list of metabolites in the GSM
>> smbl_reactions.txt - list of reactions in the GSM
>> smbl_reaction_bounds.txt - reaction bounds for reactions in the GSM
>> smbl_sij.txt - stoichiometric matrix from the GSM
>> FVA.txt - FVA results (no constraints, nutrient uptake same as in medium.txt)
>> maxbiomass_pFBA_FVA.txt - FVA results at maximal biomass value, followed by total flux limited to pFBA value at maximal biomass
>> Poplar_GPRrxns_090419.txt - tab delimited file, with reactions (column 1) and associated list of genes (column 2, comma delimited)

3) metabdatapath - stores genotype-specific metabolite levels

4) snpdatapath - stores genotype-wise sequence data


'''






####the reference genotype ID
wt_clone = '367'   ##chosen as just one gene amplicon is NA

'''load data '''

##get a list of genotypes in the study
###also copy this file in the GAMS input path 
f1 =  open(data_dir+'AllGenotypesInStudy.txt' ,'r')
temp = f1.read()
f1.close()

f1 =  open(respath_forgams+'AllGenotypesInStudy.txt' ,'w')
f1.write(temp)
f1.close()

f1 =  open(data_dir+'AllGenotypesInStudy.txt' ,'r')
data = f1.readlines()
f1.close()

cloneslist = []
for line in data:
	line = line.strip()
	cloneslist.append(line.strip('\''))
print (len(cloneslist))


###add biomass constraints to GAMS input path
f1 = open(data_dir+'Biomass_Cons_AllGenotypes.txt','r')
data = f1.read()
f1.close()

f1 = open(respath_forgams+'Biomass_Cons_AllGenotypes.txt','w')
f1.write(data)
f1.close()


f1 = open(data_dir+'Biomass_ConsList_AllGenotypes.txt','r')
data = f1.read()
f1.close()

f1 = open(respath_forgams+'Biomass_ConsList_AllGenotypes.txt','w')
f1.write(data)
f1.close()




'''
load the GSM, to define 
blocked_rxns = array, has a list of blocked reactions
ex_trans_rxns = array, has a list of exchange and transfer reactions
rxn2add_fwd = array, has a list of irreversible reactions which are in the forward direction (from pFBA + FVA results)
rxn2add_bwd = array, has a list of irreversible reactions which are in the backward direction (from pFBA + FVA results)

'''

###copy GSM files 
f1 = open(gsm_dir+'smbl_metabolites.txt','r')
data = f1.read()
f1.close()

f1 = open(respath_forgams+'metabolites.txt','w')
f1.write(data)
f1.close()

f1 = open(gsm_dir+'smbl_reaction_bounds.txt','r')
data = f1.read()
f1.close()

f1 = open(respath_forgams+'reaction_bounds.txt','w')
f1.write(data)
f1.close()

f1 = open(gsm_dir+'smbl_reactions.txt','r')
data = f1.readlines()
f1.close()

f1 = open(respath_forgams+'reactions.txt','w')
f2 = open(respath_forgams+'ExchangeRxns.txt','w')
for line in data:
	f1.write(line)
	if 'EX_' == line.strip("'")[:3]:
		f2.write(line)
f1.close()
f2.close()

####get blocked rxns , exchange, transport rxns from FVA results file
blocked_rxns = []
ex_trans_rxns = ['EX_palmitate[p]','EX_cholesterol[er]','EX_ser[p]']

####get irreversible reactions
###from the FVA runs 
rxn2add_bwd = []
rxn2add_fwd = []

f1 = open(gsm_dir + 'FVA.txt','r')
data = f1.readlines()
f1.close()
data = data[1:]
for line in data:
	line = line.strip()
	line = line.split()

	if 'TestProd_' == line[0][:9]:  ###test rxns to check for blocked biomass precursors (long-retired)
		pass
	else:
		if 'E2C_' == line[0][:4]:
			ex_trans_rxns.append(line[0])
		if 'Trans_' == line[0][:6]:
			ex_trans_rxns.append(line[0])
		if 'EX_' == line[0][:3]:
			ex_trans_rxns.append(line[0])

		###get irreversible  and blocked rxns 
		rxn = line[0]
		if rxn[-3:] == '_er':
			rxn = rxn[:-3]+'[er]'
		elif rxn[-2:] == '_c':
			rxn = rxn[:-2]+'[c]'
		elif rxn[-2:] == '_g':
			rxn = rxn[:-2]+'[g]'
		elif rxn[-2:] == '_x':
			rxn = rxn[:-2]+'[x]'
		elif rxn[-2:] == '_p':
			rxn = rxn[:-2]+'[p]'
		elif rxn[-2:] == '_m':
			rxn = rxn[:-2]+'[m]'
		else:
			#print ('weird metabolite comp format '+rxn)
			pass

		if (float(line[1]) == 0.0) and (float(line[4]) == 0.0):
			blocked_rxns.append(rxn)
		elif float(line[1]) <= 0.0 and float(line[4])< 0.0:
			rxn2add_bwd.append(rxn)

		elif float(line[1]) > 0.0 and float(line[4])>=0.0:
			rxn2add_fwd.append(rxn)

print ('# exchange + trans rxns: '+str(len(ex_trans_rxns)))
print ('# blocked: '+str(len(blocked_rxns)))
print ('# irrev fwd: ' +str(len(rxn2add_fwd)))
print ('# irrev bwd: ' + str(len(rxn2add_bwd)))

####get reaction bounds for reference flux distribution (from maxbiomass -> pFBA at max biomass -> FVA at max biomass + min network pFBA flux)
f1 = open(gsm_dir + 'maxbiomass_pFBA_FVA.txt','r')
data = f1.readlines()
f1.close()

all_rxns_bnds = []
all_rxns_sameUBLB = []
for line in data:
	line = line.strip()
	if line=='':
		pass
	elif 'TestProd_' == (line.split())[0][:9]:  ###test rxns to check for blocked biomass precursors (long-retired)
		pass
	else:
		line = line.split()
		if line[0] not in ex_trans_rxns:
			all_rxns_bnds.append(line[0])
			if float(line[1]) == float(line[4]):
				all_rxns_sameUBLB.append(line[0])

rxn_ub_dict = dict()
rxn_lb_dict = dict()

f1 = open(gsm_dir+'maxbiomass_pFBA_FVA.txt','r')
data = f1.readlines()
f1.close()

for line in data:
		line = line.strip()
		if line=='':
				pass
		else:
				rxn = line.split()[0]
				maxbnd = line.split()[1]
				minbnd = line.split()[4]
				
				rxn_lb_dict[rxn] = float(minbnd)
				rxn_ub_dict[rxn] = float(maxbnd)

print (len(list(rxn_lb_dict.keys())))

f1 = open(respath_forgams+'reaction_bounds_pFBA.txt','w')
for rxn in rxn_lb_dict:
	f1.write('LB_pfba(\''+rxn+'\') = '+str(rxn_lb_dict[rxn])+';\n')
	f1.write('UB_pfba(\''+rxn+'\') = '+str(rxn_ub_dict[rxn])+';\n')
f1.close()


####load metabolites which are expected to be in excess like water, oxygen (ref: 10.1093/bioinformatics/btw465 )
f1 = open(metabdatapath+'Metabs_Excess_list.txt','r')
metabs_excess = f1.readlines()
f1.close()

for i in range(0,len(metabs_excess)):
	metabs_excess[i] = metabs_excess[i].strip()
	metabs_excess[i] = metabs_excess[i].replace('\'','')
	metabs_excess[i] = metabs_excess[i].replace(']','')
	metabs_excess[i] = metabs_excess[i].replace('[','_')


###now we load the GSM model!
f1 = open(gsm_dir + 'smbl_sij.txt','r')
data = f1.readlines()
f1.close()

###copy in the GAMS path too
f1 = open(respath_forgams+'sij_matrix.txt','w')
for line in data:
	f1.write(line)
f1.close()

rxn2metab_dict = dict()
all_rxns = []
for i in range(0,len(data)):
	line = data[i]
	line = line.strip()
	if line =='':
		pass
	elif line[0]=='*':
		pass
	else:
		#try: 
		line_ori = line
		line = line.split()
		stoic = float(line[-1])
		line = line[0]
		line = line.split('.')
		metab = (line[0].strip()).replace('\'','')
		rxn = (line[1].strip()).replace('\'','')
		all_rxns.append(rxn)
		if stoic <0 and rxn not in rxn2add_bwd:
			if metab not in metabs_excess:  ####ignore these as pretending they are in excess
				if rxn not in rxn2metab_dict.keys():
					rxn2metab_dict[rxn] = [metab+'\t'+str(stoic)]
				else:
					rxn2metab_dict[rxn].append(metab+'\t'+str(stoic))
		if stoic >0 and rxn in rxn2add_bwd:
			if metab not in metabs_excess:
				if rxn not in rxn2metab_dict.keys():
					rxn2metab_dict[rxn] = [metab+'\t'+str(stoic)]
				else:
					rxn2metab_dict[rxn].append(metab+'\t'+str(stoic))
		#except:
		#	IndexError
		#	print (line )
		#	print ((stoic, metab, rxn ))

all_rxns = list(set(all_rxns))
print ('all gsm rxns: '+str(len(all_rxns))+'\t'+str(len(list(rxn2metab_dict.keys()))))


'''
load metabolomics data
'''


####get the metab levels for the reference genotype 
####convert metab IDs to GAMS file format 
data = np.loadtxt(metabdatapath+'MetaboliteLevels_'+wt_clone+'.tsv',dtype='str',delimiter='\t')
data = np.transpose(data)

metabID_control_ref = data[0]
metablevels_control_ref = data[1]


for i in range(0,len(metabID_control_ref)):
	metabID_control_ref[i]  = metabID_control_ref[i].strip()
	met = metabID_control_ref[i]
	met = met.replace('[','_')
	met = met.replace(']','')
	'''
	if met[:2]=='M_':
		met = met[2:]
	if met[-3:] == '_er':
		met = met[:-3]+'[er]'
	elif met[-2:] == '_c':
		met = met[:-2]+'[c]'
	elif met[-2:] == '_g':
		met = met[:-2]+'[g]'
	elif met[-2:] == '_x':
		met = met[:-2]+'[x]'
	elif met[-2:] == '_p':
		met = met[:-2]+'[p]'
	elif met[-2:] == '_m':
		met = met[:-2]+'[m]'
	else:
		print ('weird metabolite comp format '+met)
	'''
	metabID_control_ref[i] = met



##now load metabolite levels for all other genotypes, while making model constraints
f1 = open(respath_forgams+'MetabConsList_measuredMets_AllEcotypes.txt','w')
f2 = open(respath_forgams+'MetabCons_measuredMets_AllEcotypes.txt','w')

f3b = open(respath_forgams+'MetabConsList_MinMaxMets_AllEcotypes.txt','w')
f4 = open(respath_forgams+'MetabCons_MinMaxMets_AllEcotypes.txt','w')

rxns2avoid = []  ####reactions that should not be included in these constraints (useful for debugging infeasibilities)

for ecotype in cloneslist:
	
	consnum = 0
	if ecotype in cloneslist and ecotype!=wt_clone:
		f3 = open(respath+'MetabConsRatios_AllRxns_'+ecotype+'.tsv','w')
		data = np.loadtxt(metabdatapath+'MetaboliteLevels_'+ecotype+'.tsv',dtype='str',delimiter='\t')
		data = np.transpose(data)
		metabID_control = data[0]
		metablevels_control = data[1].astype(float)

		###reformat metab id
		for i in range(0,len(metabID_control)):
			metabID_control[i]  = metabID_control[i].strip()
			met = metabID_control[i]
			met = met.replace('[','_')
			met = met.replace(']','')
			metabID_control[i] = met 
			'''
			if met[:2]=='M_':
				met = met[2:]
			if met[-3:] == '_er':
				met = met[:-3]+'[er]'
			elif met[-2:] == '_c':
				met = met[:-2]+'[c]'
			elif met[-2:] == '_g':
				met = met[:-2]+'[g]'
			elif met[-2:] == '_x':
				met = met[:-2]+'[x]'
			elif met[-2:] == '_p':
				met = met[:-2]+'[p]'
			elif met[-2:] == '_m':
				met = met[:-2]+'[m]'
			else:
				print ('weird metabolite comp format '+met)

			metabID_control[i] = met 
			'''

		###divide this by control to get the final mets 
		temp_ID = []
		temp_vals = []

		for i in range(0,len(metabID_control)):
			if metabID_control[i] in metabID_control_ref:
				for j in range(0,len(metabID_control_ref)):
					if metabID_control_ref[j] == metabID_control[i]:
						val = float(metablevels_control[i])/float(metablevels_control_ref[j])
						temp_vals.append(val)
						temp_ID.append(metabID_control[i])
		###reassign 
		metabID_control = temp_ID
		metablevels_control = temp_vals
		
		rxnlist = rxn2metab_dict.keys()
		####format--> rxn: metab+'\t'+str(stoic)
		rxns_constr = []    ####list of reactions that can be constrained 	
		max_metab_ratio = -9999.0
		min_metab_ratio = 9999.0
		for rxn in rxnlist:
			temp = rxn2metab_dict[rxn]   ###metabs in current rxn 
			#print (temp)
			
			curstoic = [0]*len(temp)
			curmetabs = [0]*len(temp)
			for i in range(0,len(temp)):
				tempy = temp[i].strip()
				tempy = tempy.split('\t')
				curmetabs[i] = tempy[0]
				curstoic[i] = tempy[1]
		
			matched_control = 0
			not_matched_control = []
			
			matched_ref = 0	
			for i in range(0,len(curmetabs)):
				if curmetabs[i] in metabID_control:
					matched_control = matched_control + 1
				else:
					not_matched_control.append(curmetabs[i])
		
			if matched_control == len(curmetabs) :
				if (rxn not in ex_trans_rxns) and (rxn not in blocked_rxns):
					if rxn in rxn2add_fwd or rxn2add_bwd:

						max_flux = 0
						min_flux = 0
						cur_flux = 1.0
						
						for i in range(0,len(curmetabs)):
							cur_val_eco = 0
							
							for j in range(0,len(metabID_control)):
								if curmetabs[i] == metabID_control[j]:
									cur_val_eco = float(metablevels_control[j])
								
							cur_ratio = cur_val_eco
							if abs(cur_ratio)<1e2 and abs(cur_ratio)>1e-2:
								if cur_ratio<min_metab_ratio:
									min_metab_ratio = cur_ratio
								if cur_ratio>max_metab_ratio:
									max_metab_ratio = cur_ratio
								
							cur_flux = cur_flux*((cur_val_eco)**abs(float(curstoic[i])))
							
						#####now write this to the file 
						if abs(cur_flux)<1e2 and abs(cur_flux)>1e-2 and 'trans_' not in rxn.lower():
							cur_flux = round(float(cur_flux),5)
							rxns_constr.append(rxn)
							if rxn in all_rxns_sameUBLB:
								
								###check if omics cons can be levied
								try:
									consnum = consnum+1
									
									f2.write('Bnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =e= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+')*('+str(rxn2omics_dict[ecotype.replace('-','_')][rxn.replace('+','')])+');\n')
									f1.write('Bnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
									f3.write(rxn+'\t'+str(cur_flux)+'\t'+str(cur_flux)+'\n')

								except:
									KeyError
									consnum = consnum+1
									
									f2.write('Bnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =e= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+');\n')
									f1.write('Bnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
									f3.write(rxn+'\t'+str(cur_flux)+'\t'+str(cur_flux)+'\n')
								
							else:
								
								try:
									consnum = consnum+1
									
									f2.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') +  MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =l= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+')*('+str(rxn2omics_dict[ecotype.replace('-','_')][rxn.replace('+','')])+');\n')
									f1.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')

									
									#f2.write(rxn+'_'+ecotype+'_MinBnd..\tv(\''+rxn+'\',\''+ecotype+'\') + S_pos(\''+rxn+'\',\''+ecotype+'\') - S_neg(\''+rxn+'\',\''+ecotype+'\') =g= 0.1*LB(\''+rxn+'\',\''+ecotype+'\')*('+str(cur_flux)+');\n')
									f2.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =g= LB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+')*('+str(rxn2omics_dict[ecotype.replace('-','_')][rxn.replace('+','')])+');\n')
									f1.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
									f3.write(rxn+'\t'+str(cur_flux)+'\t'+str(cur_flux)+'\n')
								
								except:
									KeyError
									consnum = consnum+1
									
									f2.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') +  MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =l= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+');\n')
									f1.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')

									#f2.write(rxn+'_'+ecotype+'_MinBnd..\tv(\''+rxn+'\',\''+ecotype+'\') + S_pos(\''+rxn+'\',\''+ecotype+'\') - S_neg(\''+rxn+'\',\''+ecotype+'\') =g= 0.1*LB(\''+rxn+'\',\''+ecotype+'\')*('+str(cur_flux)+');\n')
									f2.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =g= LB_pfba(\''+rxn.replace('+','')+'\')*('+str(cur_flux)+');\n')
									f1.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
									f3.write(rxn+'\t'+str(cur_flux)+'\t'+str(cur_flux)+'\n')						
						#else:
						#	print 'not writing '+rxn +'\t'+str(cur_flux)

		print (ecotype+', min:'+str(min_metab_ratio)+',max:'+str(max_metab_ratio)+' rxns cons -> '+str(len(rxns_constr)))
		min_metab_ratio = round(float(min_metab_ratio),5)
		max_metab_ratio = round(float(max_metab_ratio),5)

		for rxn in rxnlist:
			if (rxn not in rxns_constr) and (rxn not in ex_trans_rxns) and (rxn not in blocked_rxns) and (rxn  in all_rxns_bnds) and (rxn.replace('+','') not in rxns2avoid) and ( (rxn in rxn2add_fwd) or (rxn in rxn2add_bwd)):
				####make it to be ive   LB_pfba<v + S_pos + S_neg <UB_pfba
				#if rxn in rxns2avoid:
				
				###multiply by curratio 
				temp = rxn2metab_dict[rxn]   ###metabs in current rxn 
				curstoic = [0]*len(temp)
				curmetabs = [0]*len(temp)
				for i in range(0,len(temp)):
					tempy = temp[i].strip()
					tempy = tempy.split('\t')
					curmetabs[i] = tempy[0]
					curstoic[i] = tempy[1]

				min_metab_ratio_currxn = ((min_metab_ratio)**abs(float(curstoic[i])))
				max_metab_ratio_currxn = ((max_metab_ratio)**abs(float(curstoic[i])))
 
				if abs(min_metab_ratio_currxn)<1e2 and abs(min_metab_ratio_currxn)>1e-2 and abs(max_metab_ratio_currxn)<1e2 and abs(max_metab_ratio_currxn)>1e-2:
					try:
						
						consnum = consnum+1
						
						#print 'should be writing to f3'
						f4.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =l= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(max_metab_ratio_currxn)+')*('+str(rxn2omics_dict[ecotype.replace('-','_')][rxn.replace('+','')])+');\n')
						f3b.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
						
						f4.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =g= LB_pfba(\''+rxn.replace('+','')+'\')*('+str(min_metab_ratio_currxn)+')*('+str(rxn2omics_dict[ecotype.replace('-','_')][rxn.replace('+','')])+');\n')
						f3b.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')

						rxns_constr.append(rxn)
					except:
						KeyError
						consnum = consnum+1
						
						#print 'should be writing to f3'
						f4.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =l= UB_pfba(\''+rxn.replace('+','')+'\')*('+str(max_metab_ratio_currxn)+');\n')
						f3b.write('MaxBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
						
						f4.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'..\tv(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + S_pos(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') - S_neg(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') + MassActSlack(\''+rxn.replace('+','')+'\',\''+ecotype.replace('-','_')+'\') =g= LB_pfba(\''+rxn.replace('+','')+'\')*('+str(min_metab_ratio_currxn)+');\n')
						f3b.write('MinBnd_'+str(consnum)+'_'+ecotype.replace('-','_')+'\n')
						rxns_constr.append(rxn)

		print (ecotype+', added: '+str(len(rxns_constr)))

f1.close()
f2.close()
f3b.close()
f4.close()

####switch lower and upper metab ratios fr bwd rxns 


f1 = open(respath_forgams+'MetabCons_MinMaxMets_AllEcotypes.txt','r')
data = f1.readlines()
f1.close()

f1 = open(respath_forgams+'MetabCons_MinMaxMets_AllEcotypes_BwdExnsflipped.txt','w')
f2 = open(respath_forgams+'MetabConsList_MinMaxMets_AllEcotypes_BwdRxnsflipped.txt','w')

count = 0
bnds_2flip_dict = dict()
rxn2eco_dict = dict() ####records all the ecos a rxn features in 
for line in data:
		#check if rxn in bwd (from pFBA)
		line = line.strip()
		if line!='':
				startind = line.find('v(')
				endind = line.find('\',')
				rxn = line[startind:endind]
				rxn = rxn.replace('v(\'','')
				if rxn in rxn_ub_dict.keys():   ####check if model has current rxn 
						if rxn_ub_dict[rxn]<0.0 and rxn_lb_dict[rxn]<0.0:
								count = count + 1
								#get bnd 
								bnd = line[line.find(')*('):]
								bnd = bnd[2:]
								#check if '*' exist
								bnd = bnd.replace(';','')

								#evaluate to take take care of bounds with omics data 
								bnd = eval(bnd)
								geno = line[line.find(',\''):line.find('\')')]
								geno = geno.replace(',','')
								geno = geno.strip('\'')
								if rxn not in rxn2eco_dict.keys():
										rxn2eco_dict[rxn] = [geno]
								else:
										rxn2eco_dict[rxn].append(geno)

								if rxn not in bnds_2flip_dict.keys():
										#bnds_2flip_dict[rxn] = [bnd]
										bnds_2flip_dict[rxn] = {}
										bnds_2flip_dict[rxn][geno] = [bnd]

								else:	
										#bnds_2flip_dict[rxn].append(bnd)
										if geno not in bnds_2flip_dict[rxn].keys():
												bnds_2flip_dict[rxn][geno] = [bnd]
										else:
												bnds_2flip_dict[rxn][geno].append(bnd)

						else:
								f1.write(line[:]+'\n')
								###gett cons name
								f2.write(line.split('..')[0]+'\n')

'''
MaxBnd_1004_BESC_401..  v('3_1_4_11_RXN[p]','BESC_401') + S_pos('3_1_4_11_RXN[p]','BESC_401') - S_neg('3_1_4_11_RXN[p]','BESC_401') + MassActSlack('3_1_4_11_RXN[p]','BESC_401') =l= UB_pfba('3_1_4_11_RXN[p]')*(1.09759);
MinBnd_1004_BESC_401..  v('3_1_4_11_RXN[p]','BESC_401') + S_pos('3_1_4_11_RXN[p]','BESC_401') - S_neg('3_1_4_11_RXN[p]','BESC_401') + MassActSlack('3_1_4_11_RXN[p]','BESC_401') =g= LB_pfba('3_1_4_11_RXN[p]')*(0.23989);

'''
consnum = 0
for rxn in bnds_2flip_dict.keys():

		genos2add = list(set(rxn2eco_dict[rxn]))
		for geno in genos2add:
				curbnds = bnds_2flip_dict[rxn][geno]
				if float(curbnds[0])>float(curbnds[1]):
						maxbnd = curbnds[0]
						minbnd = curbnds[1]
				else:
						maxbnd = curbnds[1]
						minbnd = curbnds[0]
				consnum = consnum + 1
				maxbnd = str(round(float(maxbnd),5))
				minbnd = str(round(float(minbnd),5))
				f1.write('FlpdMaxBnd_'+str(consnum)+'_'+geno+'..\tv(\''+rxn+'\',\''+geno+'\') + S_pos(\''+rxn+'\',\''+geno+'\') - S_neg(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =l= UB_pfba(\''+rxn+'\')*('+minbnd+');\n')
				f2.write('FlpdMaxBnd_'+str(consnum)+'_'+geno+'\n')
				f1.write('FlpdMinBnd_'+str(consnum)+'_'+geno+'..\tv(\''+rxn+'\',\''+geno+'\') + S_pos(\''+rxn+'\',\''+geno+'\') - S_neg(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =g= LB_pfba(\''+rxn+'\')*('+maxbnd+');\n')
				f2.write('FlpdMinBnd_'+str(consnum)+'_'+geno+'\n')

f1.close()
f2.close()





################### $$$$$$$$$$$$$$$ ###################################
###decompose stuff now!
################### $$$$$$$$$$$$$$$ ###################################

##load gene to chromosome dictionary 
gene2chrom_dict = dict()
f1 = open( snpdatapath+ 'Ptrichocarpa_444_v3.1.gene.gff3','r')
data = f1.readlines()
f1.close()

data = data[3:]
#Chr16	phytozomev12	gene	11790788	11796962	.	-	.	ID=Potri.016G113800.v3.1;Name=Potri.016G113800;ancestorIdentifier=Potri.016G113800.v3.0
#scaffold_34	phytozomev12	gene	295934	297724	.	-	.	ID=Potri.T032755.v3.1;Name=Potri.T032755
for line in data:
	line = line.strip()
	if line.split()[2] == 'gene':
		line = line.split()
		gene = line[-1]
		if ';' in gene:
			gene = gene.split(';')
			for item in gene:
				if 'ID=' in item:
					genename = item.replace('ID=','')
					genename = genename.replace('.v3.1','')
					gene2chrom_dict[genename] = line[0].replace('_','')
		else:
			genename = gene.replace('ID=','')
			genename = genename.replace('.v3.1','')
			gene2chrom_dict[genename] = line[0].replace('_','')

####run through all SNPs to keep the ones which are found in ATLEAST two different genotypes 
snp_2_eco_dict = dict()
clone_gene_snp_dict = dict()

for clone in cloneslist:
	if clone!=wt_clone:
		print ('Reading SNP lists for clone '+clone)
		gene_snp_dict = dict()
		f1 = open(snpdatapath+'SNPs_in_gene_'+clone+'.txt','r')
		data = f1.readlines()
		f1.close()

		for line in data:
			line = line.strip()
			line = line.split('\t')
			gene = line[0][:-2]
			###now get SNPs which are comma delimited
			cursnplist = line[1].split(',')
			cur_chrom = gene2chrom_dict[gene]
			for i in range(0,len(cursnplist)):
				cursnplist[i] = cur_chrom+'_'+cursnplist[i]

			if gene not in gene_snp_dict.keys():
				gene_snp_dict[gene] = cursnplist
			
			for snp in cursnplist:
				try:
					snp_2_eco_dict[snp].append(clone)
				except:
					KeyError
					snp_2_eco_dict[snp] = [clone]

		clone_gene_snp_dict[clone] = gene_snp_dict

####filter by MAF 
###MAF of 1% corresponds to 6.65 genotypes
print ('\n\nremoving SNPs that exist in just one genotype')
print (len(snp_2_eco_dict.keys()))

snps2del = []
for key in snp_2_eco_dict.keys():
	#if len(snp_2_eco_dict[key])==1:
	if len(snp_2_eco_dict[key])>0.01*len(cloneslist):   ###MAF of 1% 
		pass
	else:
		snps2del.append(key)

for snp in snps2del:
	snp_2_eco_dict.pop(snp, None)

print ('left with '+ str(len(snp_2_eco_dict.keys()) ) + ' SNPs')


###for every ecotype 		
####now decompose the S_pos(j,l) and S_neg(j,l) terms too 
####the SNPs to keep are ones found in the snp_2_eco_dict  
fcons_1 = open(respath+'DecomposedSlackIntoSNP_ConsList_AllGenotypes.txt','w')
fcons_2 = open(respath+'DecomposedSlackIntoSNP_Cons_AllGenotypes.txt','w')
fcons_3 = open(respath+'Xexists_AllGenotypes.txt','w')

##these will come in handy later, during analysis
fmat1 = open(respath+'Rxn2GeneMapping.txt','w')     ###### dimen: j x m
fmat2 = open(respath+'Gene2SNPMapping.txt','w')     #####  dimen: m x k

####using fmat3 and fmat4 in final code
fmat3 = open(respath_forgams+'Snp2EcoMapping.txt','w')	  ##### dimen: k x l 
fmat4 = open(respath_forgams+'Rxn2Snp2EcoMapping.txt','w')	  ##### dimen: j x k x l 

###define useful datastructures
snp_2_rxn_dict = dict()  ###store SNP to reaction information
allgeneslist = [] ###all genes in model

for eco in cloneslist:
	if eco != wt_clone:
		snp_genelist = list(clone_gene_snp_dict[eco].keys())  ####gives all genes with SNPs in current eco 
		print ('\n'+eco +' with '+str(len(snp_genelist))+' genes with SNPs in them (such as '+snp_genelist[0]+')')

		####this is the list of reactions to constrain in the current ecotype 
		####for this case, include ANY reaction that has atleast one gene with a SNP
		rxns_with_snp_genes = []
		f1 = open(gsm_dir+ 'Poplar_GPRrxns_090419.txt','r')  ###tab delimited file, with genes comma delimited
		data = f1.readlines()
		f1.close()
		for line in data:
			line = line.strip()
			line = line.split('\t')
			glist = line[-1].split(',')
			for gene in glist:
				if gene in snp_genelist:
					rxns_with_snp_genes.append( (line[0].replace(']','')).replace('[','_') )

		rxns_with_snp_genes = list(set(rxns_with_snp_genes))
		
		print ('# total rxns with atleast one gene with SNP: ' +str(len(rxns_with_snp_genes))+' (such as '+rxns_with_snp_genes[0]+')')

		####get genes by rxns 
		f1 = open(gsm_dir+ 'Poplar_GPRrxns_090419.txt','r')
		data = f1.readlines()
		f1.close()

		count = 0
		#for gene in snp_genelist:
		for line in data:
			line = line.strip()
			line = line.split('\t')
			rxn = line[0].strip()
			rxn = rxn.replace(']','')
			rxn = rxn.replace('[','_')
			gpr = line[1].strip()

			if (rxn in rxns_with_snp_genes) and (rxn not in blocked_rxns) and (rxn not in ex_trans_rxns):
				####get list of rxns to constrain 
				###this is different as these are not compartmentalized 
				#rxns2cons = []
				#rxns2cons = list(set(rxns2cons))
				cur_genes_added = []

				gpr = gpr.strip()

				if len(gpr.split(','))>1:
					str2write_pos = ''
					str2write_neg = ''
					gpr = gpr.split(',')
					for i in range(0,len(gpr)):
						gpr[i] = gpr[i].strip()
						temp = gpr[i]
						if (temp!='') and ('or' not in temp) and ('and' not in temp):
							
							if (temp in snp_genelist) and (temp not in cur_genes_added):
								cur_genes_added.append(temp)
								for j in range(0,len(snp_genelist)):
									if temp == snp_genelist[j]:
										cur_snplist = clone_gene_snp_dict[eco][snp_genelist[j]]
										#print cur_snplist
										#print list(set(cur_snplist))
										for k in range(0,len(cur_snplist)):
											try:
												###checking as this SNP could have been removed in MAF filtering
												exists_check = snp_2_eco_dict[cur_snplist[k]]
											#if cur_snplist[k] in snp_2_eco_dict.keys():
												snp = cur_snplist[k]
												#print 'should be writing 2 file '
												str2write_pos = str2write_pos+'X_pos(\''+snp+'\',\''+temp.replace('.','_')+'\',\''+eco.replace('-','_')+'\') + '
												str2write_neg = str2write_neg+'X_neg(\''+snp+'\',\''+temp.replace('.','_')+'\',\''+eco.replace('-','_')+'\') + '
												fcons_3.write('X_exists(\''+snp+'\',\''+temp.replace('.','_')+'\',\''+eco.replace('-','_')+'\') = 1;\n')
												
												
												fmat2.write('Gene2SNPMap(\''+temp.replace('.','_')+'\',\''+snp+'\') = 1;\n')
												fmat3.write('Snp2EcoMap(\''+snp+'\',\''+eco.replace('-','_')+'\') = 1;\n')

												#for rxn in rxns2cons:
												fmat4.write('Rxn2Snp2EcoMap(\''+rxn+'\',\''+snp+'\',\''+eco.replace('-','_')+'\') = 1;\n')
												fmat1.write('Rxn2GeneMap(\''+rxn+'\',\''+temp.replace('.','_')+'\') = 1;\n')

												###store more info in dictionaries
												if temp.replace('.','_') not in allgeneslist:
													allgeneslist.append(temp.replace('.','_'))

											except:
												KeyError

				else:
					str2write_pos = ''
					str2write_neg = ''
					gpr = gpr.strip()

					if (gpr in snp_genelist) and (gpr not in cur_genes_added):
						cur_genes_added.append(gpr)
						for j in range(0,len(snp_genelist)):
							if gpr == snp_genelist[j]:
								cur_snplist = clone_gene_snp_dict[eco][snp_genelist[j]]
								#print 'debug 2b'
								for k in range(0,len(cur_snplist)):
									try:
										exists_check = snp_2_eco_dict[cur_snplist[k]]
									#if cur_snplist[k] in snp_2_eco_dict.keys():
										snp = cur_snplist[k]
					
										str2write_pos = str2write_pos+'X_pos(\''+snp+'\',\''+gpr.replace('.','_')+'\',\''+eco.replace('-','_')+'\') + '
										str2write_neg = str2write_neg+'X_neg(\''+snp+'\',\''+gpr.replace('.','_')+'\',\''+eco.replace('-','_')+'\') + '
										fcons_3.write('X_exists(\''+snp+'\',\''+gpr.replace('.','_')+'\',\''+eco.replace('-','_')+'\') = 1;\n')
										
										fmat2.write('Gene2SNPMap(\''+gpr.replace('.','_')+'\',\''+snp+'\') = 1;\n')
										fmat3.write('Snp2EcoMap(\''+snp+'\',\''+eco.replace('-','_')+'\') = 1;\n')

										#for rxn in rxns2cons:
										fmat4.write('Rxn2Snp2EcoMap(\''+rxn+'\',\''+snp+'\',\''+eco.replace('-','_')+'\') = 1;\n')
										fmat1.write('Rxn2GeneMap(\''+rxn+'\',\''+gpr.replace('.','_')+'\') = 1;\n')

										###store more info in dictionaries
										if gpr.replace('.','_') not in allgeneslist:
											allgeneslist.append(gpr.replace('.','_'))

									except:
										KeyError

				###write stuff now
				str2write_pos = str2write_pos.strip()
				str2write_pos = str2write_pos.strip('+')
				str2write_pos = str2write_pos.strip()

				str2write_neg = str2write_neg.strip()
				str2write_neg = str2write_neg.strip('+')
				str2write_neg = str2write_neg.strip()		

				if str2write_pos!='' and str2write_neg !='':
					#if '1104' in rxn:
					#	print 'writing for '+rxn

					#for rxn in rxns2cons:
					fcons_1.write('Decompose_'+(rxn.replace('[','')).replace(']','')+'_'+eco.replace('-','_')+'_pos\n')
					fcons_2.write('Decompose_'+(rxn.replace('[','')).replace(']','')+'_'+eco.replace('-','_')+'_pos..\tS_pos(\''+rxn+'\',\''+eco.replace('-','_')+'\') =e= '+str2write_pos+';\n')
							
					fcons_1.write('Decompose_'+(rxn.replace('[','')).replace(']','')+'_'+eco.replace('-','_')+'_neg\n')
					fcons_2.write('Decompose_'+(rxn.replace('[','')).replace(']','')+'_'+eco.replace('-','_')+'_neg..\tS_neg(\''+rxn+'\',\''+eco.replace('-','_')+'\') =e= '+str2write_neg+';\n')
							
					

fcons_1.close()
fcons_2.close()				

fcons_3.close()		

fmat1.close()
fmat2.close()
fmat3.close()
fmat4.close()

###make all SNP list, eco list, and gene list 
f1 = open(respath_forgams+'SNPList_AllEcotypes.txt','w')
for key in snp_2_eco_dict.keys():
	f1.write('\''+key+'\'\n')
f1.close()

f1 = open(respath+'GeneList_AllEcotypes.txt','w')
for gene in allgeneslist:
	f1.write('\''+gene+'\'\n')
f1.close()

f1 = open(respath+'AllGenotypesInStudy.txt','w')
for eco in cloneslist:
	if eco!=wt_clone:
		f1.write('\''+eco.replace('-','_')+'\'\n')
f1.close()



####GAMS code is stuck at generating tehe objecting, so finding out the subset of rxns that are actually present in MetabCons_measuredMets_AllEcotypes.txt
#f1 = open('MetabCons_measuredMets_AllEcotypes.txt','r')
f1 = open(respath_forgams+'MetabCons_MinMaxMets_AllEcotypes.txt','r')
data = f1.readlines()
f1.close()

rxnlist = []
rxn2eeco_dict = dict()
allgenos = []
for line in data:
		line = line.strip()
		startind = line.find('v(')
		endind = line.find('\',')
		rxn = line[startind:endind]
		rxn = rxn.replace('v(\'','')
		rxnlist.append(rxn)

		geno = line[line.find(',\''):line.find('\')')]
		geno = geno.strip()
		geno = geno.strip(',')
		geno = geno.strip('\'')
		geno = geno.strip()
		if rxn not in rxn2eeco_dict.keys():
				rxn2eeco_dict[rxn] = [geno]
		else:
				rxn2eeco_dict[rxn].append(geno)
		allgenos.append(geno)

f1 = open(respath_forgams+'MetabCons_measuredMets_AllEcotypes.txt','r')
data = f1.readlines()
f1.close()

for line in data:
		line = line.strip()
		startind = line.find('v(')
		endind = line.find('\',')
		rxn = line[startind:endind]
		rxn = rxn.replace('v(\'','')
		rxnlist.append(rxn)

		geno = line[line.find(',\''):line.find('\')')]
		geno = geno.strip()
		geno = geno.strip(',')
		geno = geno.strip('\'')
		geno = geno.strip()
		if rxn not in rxn2eeco_dict.keys():
				rxn2eeco_dict[rxn] = [geno]
		else:
				rxn2eeco_dict[rxn].append(geno)
		allgenos.append(geno)


allgenos = list(set(allgenos))

'''
###check if every rxn is cons in every ecotype
for rxn in rxn2eeco_dict.keys():
		if len(list(set(rxn2eeco_dict[rxn])))!=100:
				print '\n'+ rxn +' not in '
				for geno in allgenos:
						if geno not in rxn2eeco_dict[rxn]:
								print geno
'''

rxnlist = list(set(rxnlist))


f1 = open(respath_forgams+'RxnsWithMetabCons_AllMets.txt','w')
#f1 = open('RxnsWithMetabCons_measuredMets.txt','w')
for rxn in rxnlist:
		f1.write('\''+rxn+'\'\n')

f1.close()


####do the same for Rxn2Snp2EcoMapping.txt

f1 = open(respath_forgams+'Rxn2Snp2EcoMapping.txt','r')
data = f1.readlines()
f1.close()

rxnlist = []
for line in data:
		line = line.strip()
		startind = line.find('Rxn2Snp2EcoMap(\'')
		endind = line.find('\',')
		rxn = line[startind:endind]
		rxn = rxn.replace('Rxn2Snp2EcoMap(\'','')
		rxnlist.append(rxn)

rxnlist = list(set(rxnlist))

f1 = open(respath_forgams+'RxnsWithSNPs.txt','w')
for rxn in rxnlist:
		f1.write('\''+rxn+'\'\n')
f1.close()

###then get geno-SNP mapping and make explicit v(j,l) + S_pos(j,l) - S_neg(j,l)  + MassActSlack(j,l) =l= UB_pfba(j)
###for the rest , make v(j,l) + MassActSlack(j,l) =l= UB_pfba(j)

rxn_eco_dict = dict()   #####rxn_eco_dict[rxn][eco] = metab/SNPcons/NOTHING
#v('SHIKIMATE_5_DEHYDROGENASE_RXN[p]','GW_9761')

f1 = open(respath_forgams+'MetabCons_MinMaxMets_AllEcotypes_BwdExnsflipped.txt','r')
data = f1.readlines()
f1.close()

for line in data:
        line = line.strip()
        rxn = line[line.find('v(\''):line.find('\',\'')]
        rxn = rxn.replace('v(\'','')
        geno = line[line.find('\',\''):line.find('\')')]
        geno = geno.replace('\',\'','')
        try:
                a = rxn_eco_dict[rxn]
                #try to get geno
                try:
                        b = rxn_eco_dict[rxn][geno]
                        ##exists!!! add metab
                        rxn_eco_dict[rxn][geno].append('metab')
                except:
                        KeyError
                        #add geno
                        rxn_eco_dict[rxn][geno] = ['metab']
        except:
                KeyError
                #add rxn
                rxn_eco_dict[rxn] = {}
                ###add geno
                rxn_eco_dict[rxn][geno] = ['metab']

f1 = open(respath_forgams+'MetabCons_measuredMets_AllEcotypes.txt','r')
data = f1.readlines()
f1.close()

for line in data:
        line = line.strip()
        rxn = line[line.find('v(\''):line.find('\',\'')]
        rxn = rxn.replace('v(\'','')
        geno = line[line.find('\',\''):line.find('\')')]
        geno = geno.replace('\',\'','')
        try:
                a = rxn_eco_dict[rxn]
                #try to get geno
                try:
                        b = rxn_eco_dict[rxn][geno]
                        ##exists!!! add metab
                        rxn_eco_dict[rxn][geno].append('metab')
                except:
                        KeyError
                        #add geno
                        rxn_eco_dict[rxn][geno] = ['metab']
        except:
                KeyError
                #add rxn
                rxn_eco_dict[rxn] = {}
                ###add geno
                rxn_eco_dict[rxn][geno] = ['metab']

##now read the rxn2snp2eco mapping
f1 = open(respath_forgams+'Rxn2Snp2EcoMapping.txt','r')
data = f1.readlines()
f1.close()

for line in data:
        line = line.strip()
        line = line.replace('Rxn2Snp2EcoMap(','')
        line = line.replace(') = 1;','')
        line = line.replace('\'','')
        line = line.split(',')
        rxn = line[0]
        geno = line[-1]
        try:
                a = rxn_eco_dict[rxn]
                #try to get geno
                try:
                        b = rxn_eco_dict[rxn][geno]
                        ##exists!!! add metab
                        rxn_eco_dict[rxn][geno].append('SNP')
                except:
                        KeyError
                        #add geno
                        rxn_eco_dict[rxn][geno] = ['SNP']
        except:
                KeyError
                #add rxn
                rxn_eco_dict[rxn] = {}
                ###add geno
                rxn_eco_dict[rxn][geno] = ['SNP']

f1 = open(respath_forgams+'reactions.txt','r')
allrxns = f1.readlines()
f1.close()

for i in range(0,len(allrxns)):
        allrxns[i] = allrxns[i].strip()
        allrxns[i] = allrxns[i].strip('\'')

f1 = open(respath_forgams+'AllGenotypesInStudy.txt','r')
allgenos = f1.readlines()
f1.close()

for i in range(0,len(allgenos)):
        allgenos[i] = allgenos[i].strip()
        allgenos[i] = allgenos[i].strip('\'')

f1 = open(respath_forgams+'RxnBnds_MetSNPseparated.txt','w')
f2 = open(respath_forgams+'RxnBndsList_MetSNPseparated.txt','w')
consnum = 0
for rxn in rxn_eco_dict.keys():
        if rxn.strip()!='':
                curgenos = rxn_eco_dict[rxn].keys()
                #for geno in curgenos:
                for geno in allgenos:
                        if geno in curgenos:
                                SNPcons = False
                                generalcons = False
                                if 'metab' not in rxn_eco_dict[rxn][geno]:
                                        if 'SNP' not in rxn_eco_dict[rxn][geno]:
                                                #pass
                                                consnum = consnum+1

                                                f1.write('RxnBnd1_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =l= UB_pfba(\''+rxn+'\');\n')
                                                f2.write('RxnBnd1_'+str(consnum)+'\n')
                                                f1.write('RxnBnd2_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =g= LB_pfba(\''+rxn+'\');\n')
                                                f2.write('RxnBnd2_'+str(consnum)+'\n')
                                        else:
                                                consnum = consnum+1
                                                f1.write('RxnBnd1_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + S_pos(\''+rxn+'\',\''+geno+'\') - S_neg(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =l= UB_pfba(\''+rxn+'\');\n')
                                                f2.write('RxnBnd1_'+str(consnum)+'\n')
                                                f1.write('RxnBnd2_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + S_pos(\''+rxn+'\',\''+geno+'\') - S_neg(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =g= LB_pfba(\''+rxn+'\');\n')
                                                f2.write('RxnBnd2_'+str(consnum)+'\n')

                        else:
                                #pass
                                #add general bunds
                                consnum = consnum+1
                                f1.write('RxnBnd1_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =l= UB_pfba(\''+rxn+'\');\n')
                                f2.write('RxnBnd1_'+str(consnum)+'\n')
                                f1.write('RxnBnd2_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =g= LB_pfba(\''+rxn+'\');\n')
                                f2.write('RxnBnd2_'+str(consnum)+'\n')

for rxn in allrxns:
        if rxn.strip()!='' and 'ATPM' not in rxn and 'iomass' not in rxn:
                if rxn not in rxn_eco_dict.keys():
                        for geno in allgenos:
                                #add general bunds
                                consnum = consnum+1
                                f1.write('RxnBnd1_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =l= UB_pfba(\''+rxn+'\');\n')
                                f2.write('RxnBnd1_'+str(consnum)+'\n')
                                f1.write('RxnBnd2_'+str(consnum)+'..\tv(\''+rxn+'\',\''+geno+'\') + MassActSlack(\''+rxn+'\',\''+geno+'\') =g= LB_pfba(\''+rxn+'\');\n')
                                f2.write('RxnBnd2_'+str(consnum)+'\n')

f1.close()
f2.close()
