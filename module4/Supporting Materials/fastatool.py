# Author: Garret Kern
#
# Modifications by Josh Cherry  2/16/19:
#
# 1. Fixed bug affecting cleaning of csv files.
# 2. Show user example tag field values (two for each, taken from first
#    and last sequence).
# 3. When exiting upon error, wait for user to press return before exiting.
#    Also, exit 0 for success, 1 for failure
# 4. Allow .fa in addition to .fas, and check that these are at end of fname.
# 5. Cosmetic changes: spelling, etc.

# Imports
import time
import sys
import random
import math


# Quit script after short delay
def quit_success():
    print("---Script Done---")
    time.sleep(4)
    sys.exit(0)

def quit():
    _ = input('Error in processing; press RETURN to exit: ')
    sys.exit(1)
              
# Read file for input
def inputFile(name):
    pathname = input("Filepath for " + name
                     + " (must be .fa (or .fas) or .csv format): ")
    try:
        # inputFile  = open(pathname, "r", encoding="utf8")
        inputFile  = open(pathname, "r")
    except Exception:
        print("File does not exist")
        quit()
    try:
        inputLines = inputFile.readlines()
    except Exception:
        inputFile.close()
        inputFile  = open(pathname, "r", encoding="utf8")
        inputLines = inputFile.readlines()

    inputFile.close()

    if (".csv" in pathname):
        splitChar = ","
    elif (pathname.endswith('.fa') or pathname.endswith('.fas')):
        splitChar = "|"
    else:
        print("Invalid file extension, must be \
        .csv, .fa, or .fas (the last two equivalent)")
        quit()
    return (inputLines,splitChar)


# Output to a File
def outputFile(outputLines):
    pathname = input("Filepath for output (will create a new file if none exist): ")
    outputFile = open(pathname, "w")
    for line in outputLines:
        outputFile.write(line)
    outputFile.close()


# Remove illgal fasta chars
def replaceIllegalChars(toChange):
    return toChange.rstrip().replace(":","").replace("(","_").replace(")","_").replace(" ","_").replace("\'","").replace("\"", "")


# Format input onto one line
# Can take lines from .fas or .csv file
def formatInput(inputLines, splitChar):

    # Reformat for analysis
    newLines = []
    curSequence = replaceIllegalChars(inputLines[0])
    count = 0
    for line in inputLines[1:]:
        if (">" in line):
            newLines.append(curSequence + "\n")
            curSequence = replaceIllegalChars(line)
            count = 0
        else:
            if (count == 0):
                curSequence = curSequence + splitChar
            curSequence = curSequence + line.rstrip()
            count+=1
    # Add final line, not added by for loop
    newLines.append(curSequence + "\n")
    inputLines = newLines
    return inputLines


# Format output to .fas format
def outputToFas(outputLines, splitChar):
    # Return to fasta format
    formatOutput = []
    for line in outputLines:
        tags = "|".join(line.rstrip().split(splitChar)[:-1])
        sequence = line.rstrip().split(splitChar)[-1]
        formatOutput.append(tags + "\n")
        count = 0
        chunk = ""
        for char in sequence:
            count += 1
            chunk = chunk + char
            if count == 60:
                formatOutput.append(chunk + "\n")
                chunk = ""
                count = 0
        if len(chunk) != 0:
            formatOutput.append(chunk + "\n")

    return formatOutput


class Review:

    def findMedianLength(inputLines):

        lens = []

        for line in inputLines:
            tagLen = len(line.split(Review.splitChar))
            lens.append(tagLen)

        return lens[int(len(lens)/2)]


    def findBadChars(inputLines):
        validLines = []

        for line in inputLines:
            if "\"" in line:
                print("Line: " + str(line.split(Review.splitChar)[0]) + " contains a quotation")
            elif "(" in line:
                print("Line: " + str(line.split(Review.splitChar)[0]) + " contains a paren")
            elif " " in line:
                print("Line: " + str(line.split(Review.splitChar)[0]) + " contains a space")
            elif ("," in line) and (Review.splitChar == "|"):
                print("Line: " + str(line.split(Review.splitChar)[0]) + " contains a comma")
            else:
                validLines.append(line)

        return validLines


    def findLongLines(inputLines):
        median = Review.findMedianLength(inputLines)
        validLines = []

        for line in inputLines:
            tagLen = len(line.split(Review.splitChar))
            if not(tagLen == median):
                print("Line: " + str(line.split(Review.splitChar)[0]) + " has an incorrect number of tags")
            else:
                validLines.append(line)

        return validLines


    def start():
        print("*** Starting File Reviewer ***")
        print("** Find lines with the wrong number of tags and illegal characters **")

        # Read input
        (inputLines, Review.splitChar) = inputFile("input")
        inputLines = formatInput(inputLines, Review.splitChar)

        inputLines = Review.findBadChars(inputLines)
        inputLines = Review.findLongLines(inputLines)

        # format output
        outputLines = outputToFas(inputLines, Review.splitChar)
        # Print good lines to output file
        outputFile(outputLines)

# Cleanup .fas or .csv file
class Cleanup:

    splitChar = ""

    def readCarefully(inputLines):
        threshold = float(input("Percent of total lines that must have a nucelotide in the column (enter .7 for 70%): "))
        numLines = 0
        sequenceDict = {}
        curCharDict = {}


        for line in inputLines:
            if (">" in line):
                id = line
                sequenceDict[id] = ""
                numLines += 1

        colCount = 0
        repeat = True

        while(repeat):
            # repeat as long as any lines have nucleotides in the colCount position
            repeat = False
            # count the number of valid nucleotides across all lines in a specific column
            nucCount = 0

            for line in inputLines:
                if (">" in line):
                    lineCounter = 0
                    id = line
                else:
                    lineCounter += len(line)
                    if not(colCount > lineCounter):
                        curChar = line[colCount - lineCounter]
                        if (curChar != "-"):
                            nucCount += 1
                            repeat = True
                            curCharDict[id] = curChar

            if (nucCount >= threshold * numLines):
                for line in inputLines:
                    if (">" in line):
                        id = line
                        sequenceDict[id] += curCharDict[id]
            colCount += 1

        inputLines = []

        for sequenceId in sequenceDict.keys():
            print(sequenceId + Cleanup.splitChar + sequenceDict[sequenceId])
            inputLines.append(sequenceId + Cleanup.splitChar + sequenceDict[sequenceId])

        return inputLines

    # Remove duplicates
    def removeDuplicates(inputLines):
        lineTags = [] # list of unique sequence names
        removedLines = [] # removed lines recorded in case needed for output
        newLines = [] # each unique line from inputLines

        # Show example tag field values to aid user
        def processTag(tag):
            '''
            If tag length above threshold, replace with
            truncation plus '...'
            '''
            if len(tag) > 30:
                return tag[:27] + '...'
            return tag

        print('Example tag values:')
        print()
        tags1 = inputLines[0].split(Cleanup.splitChar)
        tags2 = inputLines[-1].split(Cleanup.splitChar)
        if len(tags1) != len(tags2):
            print('Number of fields in first and last entries not equal')
            quit()
        for i in range(len(tags1)):
            print('%d  %-30s  %-30s'
                  % (i, processTag(tags1[i]), processTag(tags2[i])))
        print()

        try:
            colsComp = list(map(int, input("Enter the tag numbers for which duplicates will be compared by (starts at 0, comma-separated): ").split(",")))
        except Exception:
            print("Please type numbers separated by commas")
            quit()

        for line in inputLines:
            # Include everything but the sequence
            # Removes sequences with exactly duplicate tags
            tags = ""
            for col in colsComp:
                tags += line.rstrip().split(Cleanup.splitChar)[col]

            if tags in lineTags:
                removedLines.append(tags)
                # uncomment to print duplicates
                # print(line.rstrip().split(Cleanup.splitChar)[0])
            else:
                newLines.append(line)
                lineTags.append(tags)

        print("REMOVED DUP COUNT: " + str(len(removedLines)) )
        return newLines

    # remove rows of low length
    def removeSparseRows(inputLines):
        # Remove length
        maxLen = 0
        removedLength = []
        newLines = []

        threshold = float(input("Percent under max length for which lines will be cut (enter .75 for 75%): "))
        for line in inputLines:
            curMax = len(line.rstrip().split(Cleanup.splitChar)[-1].replace("-","").replace("n",""))
            if (curMax > maxLen):
                maxLen = curMax

        for line in inputLines:
            sequence = line.rstrip().split(Cleanup.splitChar)[-1].replace("-","").replace("n","")
            if len(sequence) < maxLen * threshold:
                removedLength.append(line)
            else:
                newLines.append(line)

        print("REMOVED LENGTH COUNT: " + str(len(removedLength)))
        return newLines

    # remove columns with too few nucleotides
    def removeSparseCols(inputLines):
        # removed columns
        threshold = float(input("Percent under max length for which columns will be cut (enter .5 for 50%): "))
        countColSize = {}
        removedColumns = []
        newLines = []
        numLines = len(inputLines)
        threshold = threshold * numLines

        # Loop through every line, count column sizes
        for line in inputLines:
            for colNum, char in enumerate(line.rstrip().split(Cleanup.splitChar)[-1]):
                    if not(colNum in countColSize.keys()):
                        countColSize[colNum] = 0
                    if not(char == "-" or char == "n"):
                        countColSize[colNum] += 1

        validCols = []
        for key in countColSize.keys():
            if not(countColSize[key] < threshold):
                validCols.append(key)

        # Loop through every line, remove columns under threshold
        for line in inputLines:
            line = line.rstrip()
            tags = Cleanup.splitChar.join(line.split(Cleanup.splitChar)[:-1])
            sequence = line.split(Cleanup.splitChar)[-1]
            validChars = [sequence[x] for x in validCols]
            newLines.append(tags + Cleanup.splitChar
                            + "".join(list(map(str, validChars))))

        removedCount = max(countColSize.keys()) - len(validCols) + 1
        print("REMOVED COLUMNS COUNT: " + str(removedCount))
        return newLines

    # Apply all cleanup functions
    def start():
        print("*** Starting Fasta Clean ***")
        print("** Remove duplicate sequences, remove sparse columns, and remove sparse sequences **")

        # Read input
        (inputLines, Cleanup.splitChar) = inputFile("input")
        inputLines = formatInput(inputLines, Cleanup.splitChar)
        print("INITIAL COUNT: " + str(len(inputLines)))

        # remove sparse cols, sparse rows, and duplicates
        outputLines = Cleanup.removeSparseCols(inputLines)
        outputLines = Cleanup.removeSparseRows(outputLines)
        outputLines = Cleanup.removeDuplicates(outputLines)

        # format output
        Cleanup.outputLines = outputToFas(outputLines, Cleanup.splitChar)

        # Print good lines to output file
        print('\a')
        outputFile(Cleanup.outputLines)


# Change .fas to .csv or .csv to .fas
class Extension:

    def toCSV(inputLines):
        outputLines = []
        splitChar = ","
        # Reformat lines
        curSequence = inputLines[0].rstrip().replace(":","").replace("(","_").replace(")","_").replace(" ","_").replace("\'","").replace("|",",")
        count = 0
        for line in inputLines[1:]:
            if (">" in line):
                outputLines.append(curSequence  + "\n")
                curSequence = line.rstrip().replace(":","").replace("(","_").replace(")","_").replace(" ","_").replace("\'","").replace("|",",")
                count = 0
            else:
                if (count == 0):
                    if (curSequence[-1] == ","):
                        curSequence = curSequence
                    else:
                        curSequence = curSequence + splitChar
                curSequence = curSequence + line.rstrip()
                count += 1
        # Add final line, not added by for loop
        outputLines.append(curSequence)
        return outputLines

    def toFas(inputLines):
        outputLines = []
        splitChar = ","
        for line in inputLines:
            tags = "|".join(line.rstrip().split(splitChar)[:-1])
            sequence = line.rstrip().split(splitChar)[-1].rstrip()
            outputLines.append(tags + "\n")
            count = 0
            chunk = ""
            for char in sequence:
                count += 1
                chunk = chunk + char
                if count == 60:
                    outputLines.append(chunk + "\n")
                    chunk = ""
                    count = 0
            if len(chunk) != 0:
                outputLines.append(chunk + "\n")
        return outputLines

    def start():
        print("*** Starting File Type Switcher ***")
        print("** Enter .csv file to change format to .fas, enter .fas file to change format to .csv **")

        (inputLines, splitChar) = inputFile("input")
        if (splitChar == "|"):
            outputLines = Extension.toCSV(inputLines)
        elif (splitChar == ","):
            outputLines = Extension.toFas(inputLines)
        else:
            print("DEBUG: impossible, file must be .fas or .csv")
            quit()

        outputFile(outputLines)


class Tag:

    def start():
        print("*** Starting Modify Tags ***")
        print("** Input a CSV, first column will be matched to a fasta sequence tag **")
        print("** Second column will be added to each matching sequence **")
        print("** Input csv should have two columns:  **")

        (inputLines, splitChar) = inputFile("sequences fasta")
        inputLines = formatInput(inputLines, Cleanup.splitChar)

        matchLines = inputFile("tags csv")[0]

        # Change tags
        matchTagNum = int(input("Which number tag is being matched from the sequences (starting at 0): "))
        matchDict = {}
        outputLines = []

        for line in matchLines:
            line = line.rstrip().split(",")
            try:
                matchDict[line[0]] = line[1]
            except Exception:
                print("Input csv must have two columns")
                quit()

        for line in inputLines:
            line = line.rstrip()
            tags = line.split(splitChar)
            toChange = str(tags[matchTagNum])
            if (toChange in matchDict.keys()):
                toAdd = matchDict[toChange]
                tags[matchTagNum] = toChange + splitChar + toAdd
                line = splitChar.join(tags)
                outputLines.append(line)
            else:
                print("Missing key for: " + tags[0] + " (key is " + toChange + ")")

        outputLines = outputToFas(outputLines, splitChar)
        outputFile(outputLines)

class Genome:

    def start():
        print("*** Starting genome data set ***")
        print("** Input any number of fasta files **")
        print("** Enter the tag number for comparison")
        print("** Output will be a csv with each tag in the inputted number spot that is in a sequence in every file ***")

        try:
            matchTag = int(input("Which tag number should be compared (starting at 0): "))
        except Exception:
            print("Tag must be an integer")
            quit()

        (inputLines, splitChar) = inputFile("input")
        inputLines = formatInput(inputLines, splitChar)
        tagList = []

        for line in inputLines:
            tag = line.rstrip().split(splitChar)[matchTag]
            tagList.append(tag)

        cont = input("Add another file (y/n): ")

        while (cont == "y"):
            (tempLines, tempSplitChar) = inputFile("input")
            tempLines = formatInput(tempLines, tempSplitChar)
            tempTagList = []
            validTags = []

            for line in tempLines:
                tag = line.rstrip().split(splitChar)[matchTag]
                tempTagList.append(tag)
            for tag in tagList:
                if tag in tempTagList and not(tag in validTags):
                    validTags.append(tag)

            tagList = validTags

            cont = input("Add another file (y/n): ")

        outputFile(",genome\n".join(tagList) + ",genome")

class Subsample:

    def start():
        print("*** Starting Sequence Subsampler Application ***")
        print("** Check readme for use instructions and examples **")

        # todo move these
        #print("|** NOMENCLATURE: tags are larger classifications, subtags are the values within the tag (Ex. Location is a tag, Canada is a subtag)")
        #print("|** For each of the following inputs the tags and subtags are inputed by the user")
        #print("|* ALL tags means every ouput sequence must have these tags (Ex. Every sequence must be from Canada )")
        #print("|* OR tags means every output sequence must be in a single OR tag group, each of which has a specified max size (Ex. 5 sequences from every year)")
        #print("|* MIN tags means if possible there will be minimum number of sequences in each MIN tag group without conflicting with ALL or OR (Ex. 1 sequence per year for each host type)")
        #print("|* NOT tags means no output sequence can have these tags (Ex. No sequences with Canada as location)")


        # Find file to be subsampled (must be a csv with tags in the first row)
        pathname = input("Filepath to be subsampled (must be .csv with a head in the first row): ")
        try:
            inputFile  = open(pathname, "r")
        except Exception:
            print("File does not exist")
            quit()

        inputLines = inputFile.readlines()
        inputFile.close()

        # Stores all tags in first line of file
        allTags = [];
        for count,tag in enumerate(inputLines[0].rstrip().split(",")):
            allTags.append((count,tag))

        # Number of each tag
        tagNums = [x[0] for x in allTags]
        tagNames = [x[1] for x in allTags]
        print("TAGS FOUND: ")
        print(*tagNames, sep="\n")

        # Remove tag line and randomize input so output will be random
        header = inputLines[0]
        inputLines = inputLines[1:]
        random.shuffle(inputLines)

        # Ask user for which tags they wish to use
        tagTypes = input("Which subsampling types would you like (of ALL,OR,MIN,NOT)(comma-separated): ").split(",")

        # Choose tags which every output sequences must have
        if ("ALL" in tagTypes):
            andTags = input("ALL tags (subtags choosen after)(comma-separated): ").split(",")
            if len([x for x in andTags if not(x in tagNames)]) > 0:
                print("One of those tags is not valid")
                quit()
        else:
            andTags = []


        andTagNums = [x[0] for x in allTags if x[1] in andTags]
        andSubTags = {} # A dictionary mapping tags to the specified subtags
        # No option for all because that is equal to leaving out the tag
        for curTag in andTags:
            if curTag == "Date" or curTag == "date":
                try:
                    startDate = [int(x) for x in input("Enter start date (year-month-day): ").split("-")]
                    endDate = [int(x) for x in input("Enter end date (year-month-day): ").split("-")]
                except Exception:
                    print("Date must be year-month-day")
                    quit()
            else:
                curInput = input("Subtags of " + curTag + " (comma-separated): ")
                andSubTags[curTag] = curInput.split(",")

        # Choose tags for which each all outputs must have one of the subtags
        # Will later choose the maximum number each group of subtags can have
        if ("OR" in tagTypes):
            orTags = input("OR tags (subtags choosen after)(comma-separated): ").split(",")
            if len([x for x in orTags if not(x in tagNames)]) > 0:
                print("One of those tags is not valid")
                quit()
        else:
            orTags = []

        orTagNums = [x[0] for x in allTags if x[1] in orTags]
        orSubTags = {}
        for tag in orTags:
            if tag == "Date" or tag == "date":
                orDateType = input("Date grouped by year (type year) or month (type month): ")
                if not(orDateType == "month" or orDateType == "year"):
                    print("Invalid date type, script will no longer function properly")
            else:
                orSubTags[tag] = input("Subtags of " + tag + " (comma-separated or type ALL): ").split(",")


        # Choose the maximum number of sequences each group can have
        # Ex: 5 tags from Canada and 5 tags from UnitedStates
        if ("OR" in tagTypes):
            maxNumSequences = input("Number per OR group: ")
            try:
                temp = int(maxNumSequences)
            except Exception:
                print("Must be an integer")
                quit()

        # Choose tags for which each all outputs must have one of the subtags
        # Will later choose the maximum number each group of subtags can have
        if ("OR" in tagTypes and "MIN" in tagTypes):
            minTags = input("MIN tags (subtags choosen after)(comma-separated): ").split(",")
        else:
            minTags = []

        minTagNums = [x[0] for x in allTags if x[1] in minTags]
        minSubTags = {}
        for tag in minTags:
            if tag == "Date" or tag == "date":
                minDateType = input("Date grouped by year (type year) or month (type month): ")
                if not(minDateType == "month" or minDateType == "year"):
                    print("Invalid date type, script will no longer function properly")
                    quit()
            else:
                minSubTags[tag] = input("Subtags of " + tag + " (comma-separated or type ALL): ").split(",")

        # Choose the minimum number of sequences each group can have
        if ("OR" in tagTypes and "MIN" in tagTypes):
            minNumSequences = input("Number per MIN group: ")
            try:
                temp = int(minNumSequences)
            except Exception:
                print("Must be an integer")
                quit()
            print("NUMBER CHOOSEN: " + minNumSequences)


        # Choose tags which every no output sequences can have
        if ("NOT" in tagTypes):
            notTags = input("NOT tags (subtags choosen after)(comma-separated): ").split(",")
            if len([x for x in notTags if not(x in tagNames)]) > 0:
                print("One of those tags is not valid")
                quit()
        else:
            notTags = []


        notTagNums = [x[0] for x in allTags if x[1] in notTags]
        notSubTags = {} # A dictionary mapping tags to the specified subtags
        # No option for all because that is equal to leaving out the tag
        for curTag in notTags:
            if curTag == "Date" or curTag =="date":
                print("No support for date not tags")
            else:
                curInput = input("Subtags of " + curTag + " (comma-separated): ")
                notSubTags[curTag] = curInput.split(",")

        # Dictionary mapping subtags to a count of how many choosen with that tag
        # Used to ensure proper number of each orTag
        orSubTagDict = {}
        minSubTagDict = {}

        # List that will contain every output line
        validLines = [];

        # Loop through every line
        for line in inputLines:
            # Strip special character such as \n from end of the line
            line = line.rstrip()
            lineTags = line.split(",")

            # Consider each line to be valid output until proven otherwise
            # Innocent until proven guilty
            valid = True

            # Key which combines different tags for groups
            # For example canada2018, unitedstates2018, canada2019
            orFullTag = ""
            # Used for the minimum groups
            minFullTag = ""
            minValid = True

            # Loop through each subtag of the line
            for tagNum, tag in enumerate(lineTags):

                # If tag was choosen as notTag check that the subTag matches
                if tagNum in notTagNums:
                    # Make sure subTag was one of input, otherwise it is not valid output
                    possibleSubTags = notSubTags[allTags[tagNum][1]]
                    if (tag in possibleSubTags):
                        valid = False

                # If tag was choosen as andTag check that the subTag matches
                if tagNum in andTagNums:

                    # If date do menial string manipulation
                    if (tagNames[tagNum] == "Date" or tagNames[tagNum] == "date"):
                        curDate = tag.split("-")
                        for i, x in enumerate(curDate):
                            if (x == ''):
                                curDate[i] = -1
                            else:
                                try:
                                    curDate[i] = int(x)
                                except Exception:
                                    print("ERROR: a entry in the date column is not yyyy-mm-dd, it is: " + x)
                                    quit()

                        correctYear = (curDate[0] >= startDate[0] and curDate[0] <= endDate[0])

                        if not(correctYear):
                            valid = False
                        elif (correctYear and not(curDate[0] > startDate[0])):
                            #if current year is the start year
                            if (curDate[1] == -1):
                                valid = False
                            elif (curDate[1] == startDate[1]):
                                if (curDate[2] == -1):
                                    valid = False
                                elif not(curDate[2] >= startDate[2]):
                                    valid = False
                            elif not(curDate[1] > startDate[1]):
                                valid = False

                        elif (correctYear and not(curDate[0] < endDate[0])):
                            #if current year is the end year
                            if (curDate[1] == -1):
                                valid = False
                            elif (curDate[1] == endDate[1]):
                                if (curDate[2] == -1):
                                    valid = False
                                elif not(curDate[2] <= endDate[2]):
                                    valid = False
                            elif not(curDate[1] < endDate[1]):
                                valid = False

                    else:
                        # Any tag but dateTag
                        # Make sure subTag was one of input, otherwise it is not valid output
                        possibleSubTags = andSubTags[allTags[tagNum][1]]
                        if not(tag in possibleSubTags):
                            valid = False

                # If inValid ignore so not added to orSubTagDict
                if not(valid):
                    break

                # If tag was choosen as orTag
                if tagNum in orTagNums:
                    if (tagNames[tagNum] == "Date" or tagNames[tagNum] == "date"):
                        curDate = tag.split("-")
                        for i, x in enumerate(curDate):
                            if (x == ''):
                                curDate[i] = -1
                            else:
                                curDate[i] = int(x)
                        dateTag = curDate[0]*100
                        if (orDateType == "month"):
                            dateTag += curDate[1]

                        orFullTag = orFullTag + str(dateTag)
                    else:
                        # If tag is not date find possible subTags
                        possibleSubTags = orSubTags[allTags[tagNum][1]]
                        if ("ALL" in possibleSubTags or tag in possibleSubTags):
                            orFullTag = orFullTag + tag
                        else:
                            valid = False

                # If tag was choosen as a minTag
                if tagNum in minTagNums:
                    if (tagNames[tagNum] == "Date" or tagNames[tagNum] == "date"):
                        curDate = tag.split("-")
                        for i, x in enumerate(curDate):
                            if (x == ''):
                                curDate[i] = -1
                            else:
                                curDate[i] = int(x)
                        dateTag = curDate[0]*100
                        if (minDateType == "month"):
                            dateTag += curDate[1]

                        minFullTag = minFullTag + str(dateTag)
                    else:
                        # If tag is not date find possible subTags
                        possibleSubTags = minSubTags[allTags[tagNum][1]]
                        if ("ALL" in possibleSubTags or tag in possibleSubTags):
                            minFullTag = minFullTag + tag
                        else:
                            minValid = False

                if valid and minValid and not(minFullTag == ""):
                    if not(minFullTag in minSubTagDict.keys()):
                        minSubTagDict[minFullTag] = []
                    if len(minSubTagDict[minFullTag]) < int(minNumSequences):
                        minSubTagDict[minFullTag].append(line) #TODO: is this necessary?
                        # instead always add this line, only remove lines when they are greater than minsubtagdict length
                    else:
                        minValid = False


            if valid and not(orFullTag == "") and not(orFullTag in orSubTagDict.keys()):
                orSubTagDict[orFullTag] = [line]
            elif valid and not(orFullTag == "") and (len(orSubTagDict[orFullTag]) < int(maxNumSequences)):
                orSubTagDict[orFullTag].append(line)
            else:
                if valid and minValid and not(minFullTag == ""):
                    # remove a line from the same orGroup as the minNumSequences
                    removed = False
                    for curToRemove in orSubTagDict[orFullTag]:
                        # Only want to remove
                        if not(curToRemove in [v for sublist in minSubTagDict.values() for v in sublist]):
                            validLines.remove(curToRemove)
                            orSubTagDict[orFullTag].remove(curToRemove)
                            orSubTagDict[orFullTag].append(line)
                            removed = True
                            break
                    if not(removed):
                        minSubTagDict[minFullTag].remove(line)
                        valid = False
                else:
                    if not(minFullTag == "" and orFullTag ==""):
                        valid = False

            if valid:
                validLines.append(line)

        # Resort lines after randomization
        validLines.sort()

        print(str(len(validLines)) + " SEQUENCES IN SUBSAMPLE ")

        # Write each sequence to the output file
        pathname = input("Filepath for output (recommended to be csv)(will create a new file if none exist): ")
        outputFile = open(pathname, "w")
        # Uncomment for output to contain header
        # outputFile.write(header)
        for line in validLines:
            outputFile.write(line + "\n")
        inputFile.close()
        outputFile.close()

        # Write meta data to a File
        outputFile = open(".metadata.csv", "w")
        outputFile.write("Key,Size,Elements \n")
        outputFile.write("--- OR GROUPS --- \n")
        for orKey in orSubTagDict.keys():
            outputFile.write(str(orKey) + "," + str(len(orSubTagDict[orKey])) + "," + str(orSubTagDict[orKey]) + "\n")
        outputFile.write("--- MIN GROUPS --- \n")
        for minKey in minSubTagDict.keys():
            outputFile.write(str(minKey) + "," + str(len(minSubTagDict[minKey])) + "," + str(minSubTagDict[minKey]) + "\n")

        outputFile.close()


def start():
    print("---Starting Fasta Modifcation Tool---\n")
    loop = "y"

    while loop == "y":
        operation = input("What operation (subsample,clean,tag,extension,review,genome): ")
        if operation == "clean":
            Cleanup.start()
        elif operation == "extension":
            Extension.start()
        elif operation == "tag":
            Tag.start()
        elif operation == "subsample":
            Subsample.start()
        elif operation == "review":
            Review.start()
        elif operation == "genome":
            Genome.start()
        else:
            print("Invalid operation")
        print("\n---Operation Done---")
        loop = input("Perform another operation (y/n): ")
        print("") # equivalent of "\n" after input
    quit_success()


# Begin script
start()
