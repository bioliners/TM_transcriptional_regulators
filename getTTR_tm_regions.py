#!/usr/bin/python
# ~ import urllib.request
import http.client
import ssl
import json
import sys
#Author: Vadim M. Gumerov

ssl._create_default_https_context = ssl._create_unverified_context

INPUT_FILE = sys.argv[1]
TM_REGIONS_OUT_FILE = "MLTF_bac_tm_regions.fa"
TM_REGIONS_COORDINATES_FILE = "MLTF_bac_tm_regions_coordinates.txt"
LOG_FILE = "MLTF_tm_regions_log.txt"

recordCounter = 0
with open(INPUT_FILE) as f1:
	for line in f1:
		splittedLine = line.split("\t")
		if len(splittedLine) <= 1:
			continue
		genomeId = splittedLine[0]
		proteinId = splittedLine[1]
		recordCounter+=1
		print (recordCounter)
		print (genomeId)
		with open(LOG_FILE, 'a') as logFile:
			logFile.write(str(recordCounter) + "\n" + genomeId + "\n")
		while True:
			try:
				conn = http.client.HTTPSConnection("mib-jouline-db.asc.ohio-state.edu")
				conn.request("GET", "/v1/genomes/" + genomeId + "/genes?search=" + proteinId + "&fields.Aseq=tmhmm2,sequence")
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
					if gene['Aseq'] == None:
						continue
					if gene['Aseq']['tmhmm2'] == None:
						continue
					if gene["Aseq"]["tmhmm2"]['tms'] == []:
						continue
					sequence = gene["Aseq"]["sequence"]
					tmmCoordsStr = ""
					for tm in gene["Aseq"]["tmhmm2"]['tms']:
						start = tm[0]
						end = tm[1]
						#start-1 becaus arrays and strings in the script are 0-based
						tmRegion = sequence[start-1:end]
						with open(TM_REGIONS_OUT_FILE, 'a') as tmFile:
							tmCoords = "tm:" + str(start) + "-" + str(end)
							tmmCoordsStr = tmmCoordsStr + tmCoords + ";"
							tmFile.write(">" + genomeId + "|" + proteinId + "|" + tmCoords + "\n" + tmRegion + "\n")
					with open(TM_REGIONS_COORDINATES_FILE, 'a') as tmCoordFile:
						tmCoordFile.write(genomeId + "\t" + proteinId + "\t" + tmmCoordsStr + "\n")
			break
		




















