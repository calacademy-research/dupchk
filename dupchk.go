package main

// 25Oct14 JBH find duplicates by looking at first 21 bases and last 21 bases of
// each read. these are store in a two UINT64 struct which is the key to a map
// of boolean values. when we see a key (the read fingerprint of 21 first and last
// bases) for the second time we add it and the reads offset to another map
// where the value is a integer array

// first we'll just loop through and see how many fingerprint dups we have

import ("os"
		"bufio"
		"compress/gzip"
		"container/heap"		
		"fmt"
		"strconv"
		"strings"
		"time"
		"log"
		"local/util"
)

const MAXLINE = 600; const CHECK_EVERY = 100000; const WAIT_AT_LEAST = time.Millisecond*250 // 07May14 change MAXLINE to 600
type ReadFingerPrint struct {
	headBases  uint64 // first 21 characters in a read stored in a uint64
	tailBases  uint64 // last 21 characters in a read stored in a uint64
}
 
var fngrDisplay func(ReadFingerPrint) string // define a function pointer
var ntVals [256]uint64 // array to hold values for all the 256 ascii chars

func main() {
	fname, chars_to_skip, show_for_grep, nTopN, topNFasta := getopts()
	if fname == "" {
		return
	}
	
	fngrDisplay = fngrStr // set function for displaying fingerprints (could be grep str also)
	if show_for_grep {
		fngrDisplay = fngrGrepStr //fngrStr // set function for displaying fingerprints (could be grep str also)
	}
	
	// read text or gzipped version of fastqc file
	isStdIn := (fname == "-") // if reading from stdin can't depend on filesize for percentage done estimate
	var scanner *bufio.Scanner
	var file *os.File
	
	if !isStdIn {
		file, _ = os.Open(fname)
		fzipped, err := gzip.NewReader(file)
		if err == nil { // it is a gzipped file
			scanner = bufio.NewScanner(fzipped)
		} else { // it is a text file
			file.Seek(0,0) // rewind to beginning since gzip test can move file offset
			scanner = bufio.NewScanner(file)
		}
		defer file.Close() // Close the file when function exits
	} else {
		scanner = bufio.NewScanner(os.Stdin)	
	}
	
	fileSize := int64(-1)
	if !isStdIn { // not stdin we can set our file size
		fi, err := file.Stat()
		if err == nil { fileSize = fi.Size() }
	}
	
	lastDisplayTime := time.Now()
	lastMsg := ""; newMsg := "" // progress msg written to stderr
	
	numRecs := 0
	numDups := 0
	initVals()
	fingerPrints := make(map[ReadFingerPrint]int)
	
	for scanner.Scan() {
		scanner.Scan() // next line is sequence
		line := scanner.Bytes() // scanner.Text() does an allocation which we don't need since each line is just used here
		
		lngth := len(line)
		encHead := ntPack3bit(line[chars_to_skip : lngth])
		encTail := ntPack3bit(line[lngth-21 : lngth])
		fngrprt := ReadFingerPrint{encHead, encTail}
		if (fingerPrints[fngrprt] > 0) {
			numDups++
		}
		fingerPrints[fngrprt] += 1
		
		numRecs++
		if ( (numRecs % CHECK_EVERY == 0) && (time.Now().Sub(lastDisplayTime) > WAIT_AT_LEAST)) {
			newMsg = util.Comma(int64(numRecs)) + " records " + pctStr(file, fileSize)
			blanks := strings.Repeat(" ", util.Max(0, len(lastMsg) - len(newMsg)))
			newMsg += blanks
			fmt.Fprint(os.Stderr, strings.Repeat("\b", len(lastMsg)) + newMsg)
			lastMsg = newMsg
			lastDisplayTime = time.Now()
		}
		
		scanner.Scan(); scanner.Scan() // throwaway qual indicator line and quality line
	}
	
	fmt.Fprintln(os.Stderr, strings.Repeat("\b", len(lastMsg)) + util.Comma(int64(numRecs)) + 
				" records total. Dups: "+util.Comma(int64(numDups))+
				"  "+util.FloatToPctStr(float64(numDups) / float64(numRecs), 2))

	if err := scanner.Err(); err != nil { log.Fatal(err) }
	
	nTopDups := showRanges(fingerPrints, nTopN)
	
	if topNFasta {
		if isStdIn {
			fmt.Fprintln(os.Stderr, "Cannot rewind stdin to write out a read for each of the TopN dups")
		} else{
			fastaTopNReads(file, chars_to_skip, nTopDups)
		}
	}
}

func initVals() {
	for ix := range ntVals {
		ntVals[ix] = 0x5 // 5 is code for non-matching char
	}
	ntVals['A'] = 0x1
	ntVals['T'] = 0x2
	ntVals['C'] = 0x3
	ntVals['G'] = 0x4
	ntVals['a'] = 0x1
	ntVals['t'] = 0x2
	ntVals['c'] = 0x3
	ntVals['g'] = 0x4
	ntVals['N'] = 0x7
	ntVals['n'] = 0x7
}

func countToRange(count int) string {
    if count < 10 {	// handle counts 1 thru 9
        return strconv.Itoa(count)
    } else if count < 20 {
		return "10-19"
    } else if count <= 50 {
		return "20-50"
    } else if count < 100 {
		return "51-99"
    } else if count < 200 {
		return "100-199"
    } else if count < 500 {
		return "200-499"
    } else if count < 1000 {
		return "500-999"
	} else if count >= 1000 {
        return "1000+"
	}
	return strconv.Itoa(count)
}

func showRanges(dups map[ReadFingerPrint]int, nTopN int) ([]util.StrInt) {
	rangeCounts := make(map[string]int)
	rangeExamples := make(map[string]ReadFingerPrint) // holds an example for each of the range types

	// minheap for TopN dups (ie those fingerprints with most dups)
	minTopN := 0;
	hp := &util.StrCtrHeap{};
	heap.Init(hp)
	for ix := 0; ix < nTopN; ix++ { // fill heap to nTopN wth zero vals, we'll replace minValue subsequently
		strctr := &util.StrInt{Key: "", Value: 0} // closure will keep this around until heap disappears
		heap.Push(hp, strctr)
	}
	
	for fngrprt, count := range dups {
		if count > 1 {
			count--
			rngStr := countToRange(count)
			rangeCounts[ rngStr ] += 1
			
			// update TopN heap if need be
			if count > minTopN { // throwaway smallest so far and replace with this new count
				heap.Pop(hp) 
				strctr := &util.StrInt{Key: fngrStr(fngrprt), Value: count} // closure will keep this around until heap disappears
				heap.Push(hp, strctr)
				minTopN = (*hp)[0].Value
			}
			
			if rangeExamples[rngStr].headBases == 0 { // we haven't yet an example for this range
				rangeExamples[rngStr] = fngrprt
			}
		}
	}
	
	fmt.Fprintln(os.Stderr, "")
	rangeList := util.SortByValue(rangeCounts)
	example := ""
	for ix := range rangeList {
		key := rangeList[ix].Key
		example = "  e.g., "
		if rangeExamples[key].headBases != 0 {
			example += fngrDisplay(rangeExamples[key])
		}
		fmt.Fprintln(os.Stderr, strconv.Itoa(rangeList[ix].Value) + " reads have " + key + " dups\t" + example)
	}

	fmt.Fprintln(os.Stderr, "\nTop " + strconv.Itoa(nTopN) + " dups")
	nTopDups := make([]util.StrInt, hp.Len())
	for ix := len(nTopDups)-1; hp.Len() > 0; ix-- { // fill backwards so we are sorted biggest to smallest (minheap not maxheap)
		tmp := heap.Pop(hp).(*util.StrInt)
		nTopDups[ix] = *tmp
	}
	for ix := 0; ix < len(nTopDups); ix++ { // now walk forwards and out th eresults
		tmp := nTopDups[ix]
		fmt.Fprintln(os.Stderr, strconv.Itoa(tmp.Value)+" dups. Match bases: "+tmp.Key)
	}
	
	return nTopDups
}

func fastaTopNReads(file *os.File, chars_to_skip int, nTopDups []util.StrInt) {
	var scanner *bufio.Scanner
	file.Seek(0,0) // rewind to beginning
	
	fzipped, err := gzip.NewReader(file)
	if err == nil { // it is a gzipped file
		scanner = bufio.NewScanner(fzipped)
	} else { // it is a text file
		file.Seek(0,0) // rewind to beginning since gzip test can move file offset
		scanner = bufio.NewScanner(file)
	}
	
	nTopFpMap := make(map[ReadFingerPrint]int)
	
	for _, strIntDup := range nTopDups { // fill fingerprint map from StrInt arrayof structs
		fp := fngrStrToFP(strIntDup.Key)
		nTopFpMap[fp] = strIntDup.Value
	}
	
	for scanner.Scan() {
		scanner.Scan() // next line is sequence
		line := scanner.Bytes() // scanner.Text() does an allocation which we don't need since each line is just used here
		
		lngth := len(line)
		encHead := ntPack3bit(line[chars_to_skip : lngth])
		encTail := ntPack3bit(line[lngth-21 : lngth])
		fngrprt := ReadFingerPrint{encHead, encTail}
		if (nTopFpMap[fngrprt] > 0) { // it is one of our top dup reads
			fmt.Println(">" + strconv.Itoa(nTopFpMap[fngrprt]) + "_dups " + fngrStr(fngrprt) + "\n" + string(line))
			delete(nTopFpMap, fngrprt)
			if len(nTopFpMap) == 0 {
				break // we've got an example for every frequent dup, so we can stop reading records in the file
			}
		}
	}
}

func pctStr(file *os.File, fileSize int64) string {
	if fileSize <= 0 {
		return ""
	}
	fpos, _ := file.Seek(0, os.SEEK_CUR) // baroque equivalent of tell() for the file offset used here
	pct := float64(fpos) / float64(fileSize)
	return util.FloatToPctStr(pct, 2)
}

func getopts() (string, int, bool, int, bool) {
    fname := ""
    skipChars := 0 // a range of gc displays to show. E.g. 1,4 will show gc content for each position and for position 1-4 totaled, 5-8 totaled, etc.
	showForGrep := false
	nTopN := 5
	topNFasta := false // if true will write out a read corresonding to each TopN dup in fasata format
	
	for ixarg := 1; ixarg < len(os.Args); ixarg++ {
		arg := os.Args[ixarg]
        if arg == "-s" { // skip chars at begin and end
            ixarg++
            if ixarg < len(os.Args) {
                skipChars, _ = util.Atoi(os.Args[ixarg])
            }
		} else if arg == "-n" { // set number of large dup fingerprints to show
            ixarg++
            if ixarg < len(os.Args) {
                nTopN, _ = util.Atoi(os.Args[ixarg])
            }

		} else if arg == "-fa" {
			topNFasta = true
		} else if arg == "-g" {
			showForGrep = true
        } else if fname == "" { // first non-flag argument is fq filename
            fname = os.Args[ixarg]
		}
	}
	
	if fname == "" { // print help msg if no filename
		fmt.Println("Usage: dupchk [[-n <NtopDups>] [fa]] <filename.fq>\n" +
					"       Checks file for duplicate reads by looking at the first 21 bases and the last 21 bases\n" +
					"       of each read. Those that match this fingerprint are considered dups for dupchk.\n" +
					"       Optional:\n" +
					"           -n <NtopDups> sets how many of the most frequent dups are shown (default: 5).\n" +
					"           -fa writes to stdout a read in fasta format for each of the top dups shown.\n" +
					"\n\n")
	}
	
	return fname, skipChars, showForGrep, nTopN, topNFasta
}

// pack upto 21 chars in ntStr into a uint64 using 3 bits per char (so 63 of the 64 bits):
// 0 is no char, 0x1 for A, 0x2 for T, 0x3 for C, 0x4 for G, 0x7 for N
// 0x5 is for any other character (such as other UIPAC ones)
func ntPack3bit(ntStr []byte) uint64 {
	var enc uint64 // upto 21 nt chars encoded in 3bits a char
	maxChrs := util.Min(21, len(ntStr))
	for ix := 0; ix < maxChrs; ix++ {
		enc = (enc << 3) | ntVals[ ntStr[ix] ]
	}
	//fmt.Printf("%s %v\n", ntStr[0:maxChrs], enc)
	return enc
}

func ntUnpack3bit(enc uint64) string {
	seq := ""
	for i := 0; i < 21; i++	{
		switch ( enc & 0x7 ) {
		case 0x1:
			seq = "A" + seq
		case 0x2:
			seq = "T" + seq
		case 0x3:
			seq = "C" + seq
		case 0x4:
			seq = "G" + seq
		case 0x7:
			seq = "N" + seq
		case 0x5: // unknown char but we'll map it to lowercase 'n'
			seq = "n" + seq
		default :
			seq = "[" + strconv.Itoa(int(enc & 0x7)) + "]" + seq
		case 0: // if bits are 000 don't output anything
		}
		
		enc >>= 3;
	}

	return seq
}

func fngrSeqs(fngr ReadFingerPrint) (headBases string, tailbases string) {
	return ntUnpack3bit(fngr.headBases), ntUnpack3bit(fngr.tailBases)
}

func fngrStr(fngr ReadFingerPrint) (bases string) {
	hd, tl := fngrSeqs(fngr)
	return hd + " " + tl
}

func fngrGrepStr(fngr ReadFingerPrint) (bases string) {
	hd, tl := fngrSeqs(fngr)
	return "^" + hd + ".*" + tl + "$"
}

func fngrStrToFP(str string) (fp ReadFingerPrint) {
	ary := strings.Split(strings.TrimSpace(str), " ")
	if len(ary) != 2 {
		return ReadFingerPrint{}
	} else {
		return ReadFingerPrint{ headBases: ntPack3bit([]byte(ary[0])), tailBases: ntPack3bit([]byte(ary[1])) }
	}
}
