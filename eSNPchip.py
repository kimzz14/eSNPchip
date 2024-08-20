###################################################################################################################
class Tools:
    def __init__(self):
        self.rcNucl_DICT = {}
        self.rcNucl_DICT['A'] = 'T'
        self.rcNucl_DICT['T'] = 'A'
        self.rcNucl_DICT['G'] = 'C'
        self.rcNucl_DICT['C'] = 'G'
        self.rcNucl_DICT['N'] = 'N'

    def reverse_complementary(self,sequence):
        result = []
        for nucl in sequence[::-1]:
            result += [self.rcNucl_DICT[nucl]]
        return ''.join(result)

tools = Tools()
###################################################################################################################
class Probe:
    def __init__(self, probeID, flank):
        self.probeID = probeID
        self.flank = flank
        self.flank_ref = None
        self.flank_alt = None
        self.flankLen = None

        self.count_ref = 0
        self.count_alt = 0

        self.ref = None
        self.alt = None

        left, ref, alt, right = flank.replace('[', '|').replace(']', '|').replace('/', '|').split('|')
        flank_ref = left + ref + right
        flank_alt = left + alt + right
        flankLen = len(flank_ref)

        self.flank_ref = flank_ref
        self.flank_alt = flank_alt
        self.flankLen = flankLen

        self.ref = ref
        self.alt = alt
    
    def str(self):
        return '\t'.join(map(str, [self.probeID, self.flank, self.flankLen, self.count_ref, self.count_alt]))
    
###################################################################################################################
class ProbeHandler:
    def __init__(self):
        self.probe_LIST = []
        self.probeID_DICT = {}
        self.flank_DICT = {}

    def read_probeFile(self, fileName):
        probe_LIST = []
        probeID_DICT = {}
        flank_DICT = {}
        fin = open(fileName)
        for line in fin:
            probeID, flank = line.rstrip('\n').split('\t')
            probe = Probe(probeID, flank)
            probe_LIST += [probe]
            probeID_DICT[probeID] = probe

            flank = probe.flank_ref
            if not flank in flank_DICT: flank_DICT[flank] = []
            flank_DICT[flank] += [(probe, 'ref')]

            flank = tools.reverse_complementary(probe.flank_ref)
            if not flank in flank_DICT: flank_DICT[flank] = []
            flank_DICT[flank] += [(probe, 'ref')]

            flank = probe.flank_alt
            if not flank in flank_DICT: flank_DICT[flank] = []
            flank_DICT[flank] += [(probe, 'alt')]

            flank = tools.reverse_complementary(probe.flank_alt)
            if not flank in flank_DICT: flank_DICT[flank] = []
            flank_DICT[flank] += [(probe, 'alt')]
        fin.close()

        self.probe_LIST = probe_LIST
        self.probeID_DICT = probeID_DICT
        self.flank_DICT = flank_DICT
    
    def make_flankFasta(self, probeFasta, flankLen):
        fout = open(probeFasta, 'w')

        for probe in self.probe_LIST:
            if probe.flankLen != flankLen: continue

            fout.write('>' + probe.probeID + '_Ref' + '\n')
            fout.write(probe.flank_ref + '\n')

            fout.write('>' + probe.probeID + '_Alt' + '\n')
            fout.write(probe.flank_alt + '\n')

        fout.close()

    def read_kmcResult(self, fileName_LIST):
        for fileName in fileName_LIST:
            fin = open(fileName)
            for line in fin:
                flank, count = line.rstrip('\n').split('\t')
                count = int(count)

                for probe, flankType in self.flank_DICT[flank]:
                    if flankType == 'ref':
                        probe.count_ref += count
                    elif flankType == 'alt':
                        probe.count_alt += count
                    else:
                        print('bug')
            fin.close()
        
    def write_result(self, outFile):
        fout = open(outFile, 'w')
        for probe in self.probe_LIST:
            fout.write('\t'.join(map(str, [probe.probeID, probe.flank, probe.flankLen, probe.count_ref, probe.count_alt])) + '\n')
        fout.close()

    def get_probe_from(self, flank):
        return self.flank_DICT[flank]

    def get_flankLen(self):
        return sorted(list(set([probe.flankLen for probe in self.probe_LIST])), reverse=True)

###################################################################################################################
import subprocess
class KMCHandler:
    def __init__(self):
        pass

    def check(self):
        command = 'kmc'
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
        stdout_LIST = process.stdout.read().decode('utf-8').split('\n')

        if len(stdout_LIST) < 5:
            return False
        else:
            return True

    def make_kmcDB(self, inFasta, outDB, dbDir):
        command = 'kmc -k71 -fm -ci0 -r {0} {1} {2}'.format(inFasta, outDB, dbDir)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
    
    def make_simple(self, inDB1, inDB2, outDB):
        command = "kmc_tools simple {0} {1} intersect {2} -ocleft".format(inDB1, inDB2, outDB)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
    
    def make_dump(self, inDB, outFile):
        command = "kmc_tools transform {0} dump {1}".format(inDB, outFile)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
    
    def rm_kmcDB(self, inDB):
        command = "rm {0}.kmc_pre".format(inDB)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
        command = "rm {0}.kmc_suf".format(inDB)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()

###################################################################################################################
from optparse import OptionParser
from joblib import Parallel, delayed
from datetime import datetime
import gzip, glob, sys

#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-f","--fastq",action = 'store',type = 'string',dest = 'FASTQ',help = "")
parser.add_option("-p","--probe",action = 'store',type = 'string',dest = 'PROBE',help = "")
parser.add_option("-o","--prefix",action = 'store',type = 'string',dest = 'PREFIX',help = "")
(opt, args) = parser.parse_args()
if opt.FASTQ == None or opt.PROBE == None or opt.PREFIX == None:
    print('Basic usage')
    print('')
    print('     python eSNPchip.py --fastq test.fastq.gz --probe 35KProbe --prefix test')
    print('')
    sys.exit()

fastqFile = opt.FASTQ
probeFile = opt.PROBE
prefix = opt.PREFIX

#fastqFile = '/archive/kimzz14/SRA_RAW/NAAS/Triticum_aestivum/Keumgang/Keumgang.02cell.hifi.fastq.gz'
#fastqFile = '../test.fastq.gz'
#prefix = 'Keumgang.02cell'

#probeFile = '35KProbe'

tmpDir = 'tmp'
lineN = 10000
batchN = 10

kmcHandler = KMCHandler()
if kmcHandler.check() == False:
    print('[{0}] kmc command not found ...'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), flush=True)
    sys.exit()
probeHandler = ProbeHandler()
probeHandler.read_probeFile(probeFile)

def split_fastq(tmpDir, prefix, lineN):
    lineN = int(lineN*4)
    fin = gzip.open(fastqFile, 'rt')
    fout = None
    for lineIDX, line in enumerate(fin):
        if lineIDX%4 == 0:
            if lineIDX%lineN == 0:
                if fout != None: 
                    fout.close()
                fileIDX = int(lineIDX/lineN)
                fileName = '{0}/{1}-{2:05}.fa'.format(tmpDir, prefix, fileIDX)
                fout = open(fileName, 'w')
        if lineIDX%4 == 1:
            fout.write('>test' + '\n')
            fout.write(line)
    fout.close()

fileN = len(glob.glob('{0}/{1}-?????.fa'.format(tmpDir, prefix)))
if fileN == 0:
    print('[{0}] split_fastq running ...'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), flush=True)
    split_fastq(tmpDir, prefix, lineN)
else:
    print('[{0}] split_fastq skip '.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), flush=True)

def run_flankLen(flankLen):
    fileN = len(glob.glob('{0}/{1}-?????.fa'.format(tmpDir, prefix)))
    probeFasta  = '{0}/{1}-probe-K{2:03}.fa'.format(tmpDir, prefix, flankLen)
    probeDB     = '{0}/{1}-probe-K{2:03}'.format(tmpDir, prefix, flankLen)

    probeHandler.make_flankFasta(probeFasta, flankLen)
    kmcHandler.make_kmcDB(probeFasta, probeDB, tmpDir)

    def run_batch_kmc(tmpDir, prefix, flankLen, fileN, batchN):
        def run_single_kmc(fileIDX):
            probeDB  = '{0}/{1}-probe-K{2:03}'.format(tmpDir, prefix, flankLen)
            inFasta  = '{0}/{1}-fastq-{3:05}.fa'.format(tmpDir, prefix, flankLen, fileIDX)
            fastqDB  = '{0}/{1}-fastq-K{2:03}-{3:05}'.format(tmpDir, prefix, flankLen, fileIDX)
            commomDB = '{0}/{1}-fastq-K{2:03}-{3:05}-commom'.format(tmpDir, prefix, flankLen, fileIDX)
            outFile  = '{0}/{1}-fastq-K{2:03}-{3:05}-commom.txt'.format(tmpDir, prefix, flankLen, fileIDX)

            kmcHandler.make_kmcDB(inFasta, fastqDB, tmpDir)
            kmcHandler.make_simple(fastqDB, probeDB, commomDB)
            kmcHandler.make_dump(commomDB, outFile)
            kmcHandler.rm_kmcDB(fastqDB)
            kmcHandler.rm_kmcDB(commomDB)

        Parallel(n_jobs=batchN)(delayed(run_single_kmc)(i) for i in range(0, fileN))

    run_batch_kmc(tmpDir, prefix, flankLen, fileN, batchN)
    kmcHandler.rm_kmcDB(probeDB)

for flankLen in probeHandler.get_flankLen()[0:1]:
    resultN = len(glob.glob('{0}/{1}-K{2:03}-?????-commom.txt'.format(tmpDir, prefix, flankLen)))
    if resultN == 0:
        print('[{0}] run_flankLen({1}) running ...'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), flankLen), flush=True)
        run_flankLen(flankLen)
    else:
        print('[{0}] run_flankLen({1}) skip ...'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), flankLen), flush=True)

print('[{0} done ...'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), flush=True)
fileName_LIST = glob.glob('{0}/{1}-K???-?????-commom.txt'.format(tmpDir, prefix))
probeHandler.read_kmcResult(fileName_LIST)
probeHandler.write_result(prefix + '.txt')