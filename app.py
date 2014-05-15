import csv
import sys
import hgvs
import hgvs.utils
#from hgvs import utils
#import pygr
from pygr.seqdb import SequenceFileDB
from flask import Flask,jsonify, request


application = app = Flask(__name__)

# genome = SequenceFileDB('/local/resources/hg18/karyotypic/hg18.fa')
genome = SequenceFileDB('/Users/afrieden/test/hg/hg18.fa')



geneTranscripts = {'CFTR':'NM_000492.3','PCDH15':'NM_033056.3','ABCC8':'NM_000352.3','ASPA':'NM_000049.2','BCKDHA':'NM_000709.3','BCKDHB':'NM_000056.3',
'BLM':'NM_000057.2','CLRN1':'NM_174878.2','DLD':'NM_000108.3','FANCC':'NM_000136.2','G6PC':'NM_000151.3','HEXA':'NM_000520.4',
							'IKBKAP':'NM_003640.3',
							'MCOLN1':'NM_020533.2',
							'SMPD1':'NM_000543.4',
							'TMEM216':'NM_001173990.2',
							'FKTN':'NM_006731.2',
							'GBA':'NM_001005741.2',
							'NEB':'NM_001164507.1'}

def getGeneFromPosition(_position):
	#check which gene it is in
	geneExportPath = '/Users/afrieden/Projects/hgvs/hgvs/data/geneExport.csv'
	reader = csv.reader(open(geneExportPath), delimiter=',')
	next(reader, None)
	contents = [line for line in reader] 
	position = int(_position)
	for row in contents:
		if(position > int(row[1])):
			if(position < int(row[2])):
				return row[0]
	_max = sys.maxint
	closestGene = ''
	for row in contents:
		diff_start = abs(int(row[1])-position)
		diff_end = abs(int(row[2]) - position)
		diff = min(diff_start,diff_end)
		if(diff < _max):
			_max = diff
			closestGene = row[0]
	return closestGene

#infile = open('/sandbox/afrieden/Projects/hgvs/hgvs/data/gsg-transcript-03-06-2014.txt')
with open('/Users/afrieden/Projects/hgvs/hgvs/data/gsg-transcript-03-06-2014.txt') as infile:
	transcripts = hgvs.utils.read_transcripts(infile)
#transcripts = read_transcripts(infile)

def get_transcript(name):
    return transcripts.get(name)




@app.route('/')
def index():
    return jsonify( {
    	'first':'http://localhost:5000/convert/cdnaToVcf?gene=CFTR&cdna=c.1521_1523delCTT',
    	'second':'http://localhost:5000/convert/vcfToCdna?chrom=chr7&pos=116986880&ref=ATCT&alt=A'
    	})



@app.route('/convert/vcfToCdna')
def vcfToCdna():
	_chrom = str(request.args.get('chrom'))
	_pos = int(request.args.get('pos'))
	_ref = str(request.args.get('ref'))
	_alt = str(request.args.get('alt'))
	
	chrom, offset, ref, alt = (_chrom, _pos, _ref, _alt)
	gene = str(getGeneFromPosition(_pos))
	transcript = get_transcript(geneTranscripts[gene])
	hgvs_name = hgvs.format_hgvs_name(
    	chrom, offset, ref, alt, genome, transcript)
	cdnaName = hgvs_name.split(':')[1]
	fullName = gene + ':' + cdnaName
	return jsonify( {'name':fullName})



@app.route('/convert/cdnaToVcf')
def cdnaToVcf():
	gene = request.args.get('gene')
	cdna_name = request.args.get('cdna')
	transcript = geneTranscripts[gene]
	transcriptName = str(transcript + ':' + cdna_name)
	chrom, offset, ref, alt = hgvs.parse_hgvs_name(
		transcriptName, genome, get_transcript=get_transcript)

	return jsonify({ 	
		'foundTranscript': transcript,
		'geneInput':gene,
		'cdna_name':cdna_name,
		'transcriptName':transcriptName,
		'chrom':chrom,
		'position':offset,
		'ref':ref,
		'alt':alt

	})





if __name__ == '__main__':
    app.run(debug = True)
