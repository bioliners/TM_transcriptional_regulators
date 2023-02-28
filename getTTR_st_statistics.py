#!/usr/bin/python
import http.client
import ssl
import json
import sys
#Author: Vadim M. Gumerov

ssl._create_default_https_context = ssl._create_unverified_context

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]
LOG_FILE = "ST_data_log.txt"

with open(LOG_FILE, 'w') as logFile:
	logFile.write("Log of " + sys.argv[0] + " running.\n")
with open(OUTPUT_FILE, 'w') as output_file:
	output_file.write( "\t".join('GenomeId, NumOfMLTFs, Total, OCP_total, TCP_total, TCP_HK, TCP_HHK, TCP_RR, TCP_HRR, TCP_Other, ChemSys, ECF, Other, \n'.split(",")) )

recordCounter = 0
with open(INPUT_FILE) as inputFile:
	for line in inputFile:	
		splittedLine = line.split("\t")
		numOfMLTFs = splittedLine[0].strip()
		genomeId = splittedLine[1].strip()

		recordCounter+=1
		with open(LOG_FILE, 'a') as logFile:
			logFile.write(str(recordCounter) + "\n" + genomeId + "\n")
			
		while True:
			try:
				conn = http.client.HTTPSConnection("api.mistdb.com")
				conn.request("GET", "/v1/genomes/" + genomeId + "/stp-matrix?per_page=1")
				res = conn.getresponse()
				html = res.read().decode("utf-8")
				if html == '[]':
					break
				data = json.loads(html)					
				if "name" in data:	#404 NotFoundError
					break
			#except json.decoder.JSONDecodeError:	#504 Gateway timeout
			except ValueError:	#504 Gateway timeout;
				continue
			else:
				total = OCP_total = TCP_total = TCP_HK = TCP_HHK = TCP_RR = TCP_HRR = TCP_Other = ChemSys = ECF = Other = "NaN"
				if "counts" in data:
					if "$total" in data["counts"]:
						total = data["counts"]["$total"]
					if "ocp" in data["counts"]: 	
						OCP_total = data["counts"]["ocp"]["$total"]
					if "tcp" in data["counts"]: 
						TCP_total = data["counts"]["tcp"]["$total"]
						if "hk" in data["counts"]["tcp"]:
							TCP_HK = data["counts"]["tcp"]["hk"]["$total"]
						if "hhk" in data["counts"]["tcp"]:	
							TCP_HHK = data["counts"]["tcp"]["hhk"]["$total"]
						if "rr" in data["counts"]["tcp"]:	
							TCP_RR = data["counts"]["tcp"]["rr"]["$total"]
						if "hrr" in data["counts"]["tcp"]:	
							TCP_HRR = data["counts"]["tcp"]["hrr"]["$total"]
						if "other" in data["counts"]["tcp"]:		
							TCP_Other = data["counts"]["tcp"]["other"]["$total"]
					if "chemotaxis" in data["counts"]:
						ChemSys = data["counts"]["chemotaxis"]["$total"]
					if "ecf" in data["counts"]:
						ECF = data["counts"]["ecf"]["$total"]
					if "other" in data["counts"]:
						Other = data["counts"]["other"]["$total"]

				with open(OUTPUT_FILE, 'a') as output_file:
					output_file.write( "\t".join(map(str, ([genomeId, numOfMLTFs, total, OCP_total, TCP_total, TCP_HK, TCP_HHK, TCP_RR, TCP_HRR, TCP_Other, ChemSys, ECF, Other]))) + "\n")						
			break		
		



















