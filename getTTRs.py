#!/usr/bin/python
# ~ import urllib.request
import http.client
import ssl
import json
#Authors: Vadim M. Gumerov, Jiawei Xing

ssl._create_default_https_context = ssl._create_unverified_context

dna_binding_domains = 'Arc	MetJ	Omega_Repress	PSK_trans_fac	RHH_1	RHH_3	RHH_4	RHH_7	AphA_like	Arg_repressor	BetR	Cro	Crp	FaeA	Fe_dep_repress	FeoC	FUR	GcrA	GerE	GntR	HrcA_DNA-bdg	HSF_DNA-bind	HTH_1	HTH_10	HTH_11	HTH_12	HTH_13	HTH_15	HTH_17	HTH_18	HTH_19	HTH_20	HTH_21	HTH_22	HTH_23	HTH_24	HTH_25	HTH_26	HTH_27	HTH_28	HTH_29	HTH_3	HTH_30	HTH_31	HTH_32	HTH_33	HTH_34	HTH_35	HTH_36	HTH_37	HTH_38	HTH_39	HTH_40	HTH_41	HTH_43	HTH_45	HTH_46	HTH_47	HTH_5	HTH_6	HTH_7	HTH_8	HTH_AraC	HTH_AsnC-type	HTH_CodY	HTH_Crp_2	HTH_DeoR	HTH_IclR	HTH_Mga	HTH_WhiA	HxlR	KORA	KorB	LacI	LexA_DNA_bind	MarR	MarR_2	MerR	MerR-DNA-bind	MerR_1	MerR_2	Mga	Mor	PaaX	PadR	Pencillinase_R	Phage_CI_repr	PuR_N	Put_DNA-bind_N	RepL	Rrf2	SgrR_N	Sigma70_r2	Sigma70_r3	Sigma70_r4	Sigma70_r4_2	TetR_N	Trans_reg_C	TrmB	Trp_repressor	UPF0122	Vir_act_alpha_C	ArsD	ComK	CtsR	LytTR	ROS_MUCR	DNA_PPF'.split()

representative_genomes = []
#bac120_metadata_r202_rep.tsv
with open("ar122_metadata_r202_rep.tsv") as f1:
	for line in f1:
		if "representative genome" in line:
			if line.split()[0][5] == 'F':
				representative_genomes.append(line.split()[0][3:])
print(len(representative_genomes))


def removeOverlapps(domainsSorted):
	tolerance = 10
	pfamFinal = []
	pfam1 = domainsSorted[0]
	significantPfam = pfam1
	overlapLength = None
	lastAdded = pfam1
	pfamFinal.append(pfam1)
	for pfam2 in domainsSorted[1:]:
		if pfam1["ali_to"] > pfam2["ali_from"]:
			overlapLength = pfam1["ali_to"] - pfam2["ali_from"]
			if overlapLength > tolerance:
				significantPfam = compareEvalues(pfam1, pfam2)
				if lastAdded == pfam1 and lastAdded != significantPfam:
					pfamFinal.remove(lastAdded)
					pfamFinal.append(significantPfam)
				lastAdded = significantPfam
			else:
				pfamFinal.append(pfam2)
				lastAdded = pfam2
				significantPfam = pfam2
			pfam1 = significantPfam
		else:
			pfamFinal.append(pfam2)
			lastAdded = pfam2
			significantPfam = pfam2
			pfam1 = significantPfam
	return pfamFinal


def compareEvalues(pfam1, pfam2):
	if "ieValue" in pfam1:
		eval1 = pfam1["i_evalue"]
		eval2 = pfam2["i_evalue"]
	else:
		eval1 = pfam1["c_evalue"]
		eval2 = pfam2["c_evalue"]
	significantPfam = None
	if eval1 < eval2:
		significantPfam = pfam1
	elif eval1 > eval2:
		significantPfam = pfam2
	elif eval1 == eval2:
		if (pfam1["ali_to"] - pfam1["ali_from"]) >= (pfam2["ali_to"] - pfam2["ali_from"]):
			significantPfam = pfam1
		else:
			significantPfam = pfam2
	return significantPfam


for genome in representative_genomes:
	for i in range(1, 101):
		while True:
			try:
				conn = http.client.HTTPSConnection("mib-jouline-db.asc.ohio-state.edu")
				conn.request("GET", "/v1/genomes/" + genome + "/genes?fields.Aseq=pfam31,sequence,tmhmm2&page=" + str(i) + "&per_page=100")
				res = conn.getresponse()
				html = res.read().decode("utf-8")
				# ~ url = 'https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/' + genome + '/genes?fields.Aseq=pfam31,sequence,tmhmm2&page=' + str(i) + '&per_page=100'
				# ~ html = urllib.request.urlopen(url).read().decode("utf-8")
				if html == '[]':
					break
				with open('temp.json', 'w') as f2:
					f2.write(html)
				with open('temp.json') as f3:
					data = json.load(f3)
				if "name" in data:	#404 NotFoundError
					break
			#except json.decoder.JSONDecodeError:	#504 Gateway timeout
			except ValueError:	#504 Gateway timeout;
				continue
			else:
				for gene in data:
					print(gene['version'])
					if gene['Aseq'] == None:
						continue
					if gene['Aseq']['tmhmm2'] == None:
						continue
					if gene["Aseq"]["tmhmm2"]['tms'] == []:
						continue
					if gene['Aseq']['pfam31'] == []:
						continue
						
					domainsSorted = sorted(gene["Aseq"]['pfam31'], key=lambda x: x["ali_from"], reverse=False)
					domainsWithNoOverlaps = removeOverlapps(domainsSorted)

					for domain in domainsWithNoOverlaps:
						if domain['name'] in dna_binding_domains:
							with open("MLTF_archaea.fa", 'a') as f5:
								f5.write(">" + genome + "\t" + gene['version'] + '\t' + str(gene['old_locus']) + "\n" + gene["Aseq"]['sequence'] + "\n")
							with open("MLTF_archaea.tsv", 'a') as f4:
								f4.write(genome + "\t" + gene['version'] + '\t' + str(gene['old_locus']) + "\t" + str(len(gene["Aseq"]["tmhmm2"]['tms'])) + "\t")
								
								domain_names = []
								for domain1 in domainsWithNoOverlaps[:-1]:
									domain_names.append(domain1['name'])
									f4.write(domain1['name'] + ",")
								domain_names.append(domainsWithNoOverlaps[-1]['name'])
								f4.write(domainsWithNoOverlaps[-1]['name'] + '\t')
								
								domain_names.sort()
								for domain2 in domain_names[:-1]:
									f4.write(domain2 + ',')
								f4.write(domain_names[-1] + '\t' + str(len(gene['Aseq']['sequence'])) + '\n')
							break
			break

