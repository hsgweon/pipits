import sys
import uc

FileName = sys.argv[1]

def OnRec():
	global OTUs, Samples, OTUTable
	if uc.Type != 'H':
		return

	OTUId = uc.TargetLabel
	if OTUId not in OTUIds:
		OTUIds.append(OTUId)
		OTUTable[OTUId] = {}

	SampleId = uc.QueryLabel.split("_")[0]

	if SampleId not in SampleIds:
		SampleIds.append(SampleId)

	if "size=" in uc.QueryLabel:
		N = int(uc.QueryLabel.split("size=")[1].split(";")[0])
	else:
		N = 1

	try:
		OTUTable[OTUId][SampleId] += N
	except:
		OTUTable[OTUId][SampleId] = N

OTUIds = []
SampleIds = []
OTUTable = {}

uc.ReadRecs(FileName, OnRec)

s = "#OTU_ID"
for SampleId in SampleIds:
	s += "\t" + SampleId
print s

for OTUId in OTUIds:
	s = OTUId
	for SampleId in SampleIds:
		try:
			n = OTUTable[OTUId][SampleId]
		except:
			n = 0
		s += "\t" + str(n)
	print s
