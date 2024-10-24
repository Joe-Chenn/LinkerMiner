from MSTool import toolCountCharInString, toolGetWord
from MSSysterm import MARK_LABEL_INFO, VALUE_ILLEGAL
from MSLogging import logToUser, logGetWarning, logGetError, INFO_TO_USER_FunctionComposition

class CFunctionComposition:

	def __init__(self, inputDP):
		self.dp = inputDP

	def __captainGetNumberElement8Mod(self, inputMod, inputE, inputDictModCom):

		result = 0
		comMod = inputDictModCom[inputMod]
		iLastRightBracket = 0  # java可以写为-1，python不行
		iELE = VALUE_ILLEGAL

		# find ele
		for i in range(len(comMod)):

			if comMod[i] == '(':

				strELE = comMod[iLastRightBracket:i]

				if strELE == inputE:
					iELE = i - 1
					break

			if comMod[i] == ')':
				iLastRightBracket = i + 1

		# get number
		if iELE == VALUE_ILLEGAL:

			return result  # 此时result为0

		else:

			p_bracket_left = 0
			p_bracket_right = 0

			for j in range(iELE + 1, len(comMod)):

				if comMod[j] == '(':
					p_bracket_left = j

				if comMod[j] == ')':
					p_bracket_right = j
					break

			result = int(comMod[p_bracket_left + 1:p_bracket_right])

		return result

	def __captainGetNumberElement8Glycan(self, inputGlycan, inputE, inputDictGlycanCom):  # TODO 这个没改好

		result = 0
		comGlycan = inputDictGlycanCom[inputGlycan]
		iLastRightBracket = 0  # java可以写为-1，python不行
		iELE = VALUE_ILLEGAL

		# find ele
		for i in range(len(comGlycan)):

			if comGlycan[i] == '(':

				strELE = comGlycan[iLastRightBracket:i]

				if strELE == inputE:
					iELE = i - 1
					break

			if comGlycan[i] == ')':
				iLastRightBracket = i + 1

		# get number
		if iELE == VALUE_ILLEGAL:

			return result  # 此时result为0

		else:

			p_bracket_left = 0
			p_bracket_right = 0

			for j in range(iELE + 1, len(comGlycan)):

				if comGlycan[j] == '(':
					p_bracket_left = j

				if comGlycan[j] == ')':
					p_bracket_right = j
					break

			result = int(comGlycan[p_bracket_left + 1:p_bracket_right])

		return result

	def __captainGetNumberElement8Link(self, inputLink, inputE, inputDictLinkCom):

		# init
		result = 0
		comLink = inputDictLinkCom[inputLink]
		iLastRightBracket = 0  # java可以写为-1，python不行
		iELE = VALUE_ILLEGAL

		# find ele
		for i in range(len(comLink)):

			if comLink[i] == '(':

				strELE = comLink[iLastRightBracket:i]

				if strELE == inputE:
					iELE = i - 1
					break

			if comLink[i] == ')':
				iLastRightBracket = i + 1

		# get number
		if iELE == VALUE_ILLEGAL:

			return result  # 此时result为0

		else:

			p_bracket_left = 0;
			p_bracket_right = 0

			for j in range(iELE + 1, len(comLink)):

				if comLink[j] == '(':
					p_bracket_left = j

				if comLink[j] == ')':
					p_bracket_right = j
					break

			result = int(comLink[p_bracket_left + 1:p_bracket_right])

		return result

	def __captainGetNumberElement8AA(self, inputAA, inputE, inputDictAACom):

		# init
		result = 0
		comAA = inputDictAACom[inputAA]
		iLastRightBracket = 0  # java可以写为-1，python不行
		iELE = VALUE_ILLEGAL

		# find ele
		for i in range(len(comAA)):

			if comAA[i] == '(':

				strELE = comAA[iLastRightBracket:i]

				if strELE == inputE:

					iELE = i - 1
					break

			if comAA[i] == ')':

				iLastRightBracket = i + 1

		# get number
		if iELE == VALUE_ILLEGAL:

			return result  # 此时result为0

		else:

			p_bracket_left = 0;
			p_bracket_right = 0

			for j in range(iELE+1, len(comAA)):

				if comAA[j] == '(':
					p_bracket_left = j

				if comAA[j] == ')':
					p_bracket_right = j
					break

			result = int(comAA[p_bracket_left+1:p_bracket_right])

		return result

	def parseLabelInfo(self, inputStrLabelInfo, inputSeq, inputMod, inputGLC, inputLIK, inputINI):

		# NONE 或者 AA:R:N:15N&AA:R:C:13C&AA:K:C:13C&AA:K:N:15N

		result = ""

		if len(inputStrLabelInfo) == 4 and inputStrLabelInfo.upper() == MARK_LABEL_INFO[0]:

			return result

		else:

			if inputStrLabelInfo[-1] == '&':
				pass
			else:
				inputStrLabelInfo = inputStrLabelInfo + '&'

			nInfo = toolCountCharInString(inputStrLabelInfo, '&')

			for iInfo in range(nInfo):

				subStrLabelInfo = toolGetWord(inputStrLabelInfo, iInfo, '&')

				if len(subStrLabelInfo) > 3 and subStrLabelInfo[0:3].upper() == MARK_LABEL_INFO[1]:  # AA:

					tmpSeq = inputSeq + '?'  # +是专门为18O这种准备的

					markAA = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					for tmpAA in tmpSeq:

						if tmpAA == markAA or markAA == '*':

							nELE = self.__captainGetNumberElement8AA(tmpAA, eleLight, self.dp.myINI.DICT1_AA_COM)
							result = result + eleLight
							result = result + '(-'
							result = result + str(nELE)
							result = result + ')'
							result = result + eleHeavy
							result = result + '('
							result = result + str(nELE)
							result = result + ')'

				elif len(inputMod) > 0 and subStrLabelInfo[0:4].upper() == MARK_LABEL_INFO[2]:  # MOD:
				# elif subStrLabelInfo[0:4].upper() == MARK_LABEL_INFO[2]:  # MOD:

					markMod = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					nMod = toolCountCharInString(inputMod, ';')

					for iMod in range(nMod):

						tmpMod = toolGetWord(toolGetWord(inputMod, iMod, ';'), 1, ',')

						if tmpMod == markMod:
							nELE = self.__captainGetNumberElement8Mod(markMod, eleLight, self.dp.myINI.DICT2_MOD_COM)
							result = result + eleLight
							result = result + '(-'
							result = result + str(nELE)
							result = result + ')'
							result = result + eleHeavy
							result = result + '('
							result = result + str(nELE)
							result = result + ')'

				elif len(inputLIK) > 0 and subStrLabelInfo[0:5].upper() == MARK_LABEL_INFO[5]:

					markLink = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					nELE = self.__captainGetNumberElement8Link(markLink, eleLight, self.dp.myINI.DICT4_LINKER_COM)
					result = result + eleLight
					result = result + '(-'
					result = result + str(nELE)
					result = result + ')'
					result = result + eleHeavy
					result = result + '('
					result = result + str(nELE)
					result = result + ')'

				elif len(inputLIK) > 0 and toolGetWord(subStrLabelInfo, 0, ':').upper() == MARK_LABEL_INFO[3]:  # 这是专门为交联4种15N情况设计的

					markAA = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					subSeq = toolGetWord(inputSeq, 0, '-')  # 和镇霖确认过，第一个是大质量的

					for tmpAA in subSeq:

						if tmpAA == markAA or markAA == '*':
							nELE = self.__captainGetNumberElement8AA(tmpAA, eleLight, self.dp.myINI.DICT1_AA_COM)
							result = result + eleLight
							result = result + '(-'
							result = result + str(nELE)
							result = result + ')'
							result = result + eleHeavy
							result = result + '('
							result = result + str(nELE)
							result = result + ')'

				elif len(inputLIK) > 0 and toolGetWord(subStrLabelInfo, 0, ':').upper() == MARK_LABEL_INFO[
								4]:  # 这是专门为交联4种15N情况设计的

					markAA = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					subSeq = toolGetWord(inputSeq, 1, '-')  # 和镇霖确认过，第二个是小质量的

					for tmpAA in subSeq:

						if tmpAA == markAA or markAA == '*':
							nELE = self.__captainGetNumberElement8AA(tmpAA, eleLight,
																	 self.dp.myINI.DICT1_AA_COM)
							result = result + eleLight
							result = result + '(-'
							result = result + str(nELE)
							result = result + ')'
							result = result + eleHeavy
							result = result + '('
							result = result + str(nELE)
							result = result + ')'

				elif len(inputGLC) > 0 and toolGetWord(subStrLabelInfo, 0, ':').upper() == MARK_LABEL_INFO[6]:  # GLYCAN:

					markGlycan = toolGetWord(subStrLabelInfo, 1, ':')
					eleLight = toolGetWord(subStrLabelInfo, 2, ':')
					eleHeavy = toolGetWord(subStrLabelInfo, 3, ':')

					nGlycan = toolCountCharInString(inputGLC, ';')

					for iGlycan in range(nGlycan):

						tmpGlycan = toolGetWord(toolGetWord(inputGLC, iGlycan, ';'), 1, ',')

						if tmpGlycan == markGlycan or markGlycan == '*':
							nELE = self.__captainGetNumberElement8Glycan(tmpGlycan, eleLight, self.dp.myINI.DICT3_GLYCO_COM)
							result = result + eleLight
							result = result + '(-'
							result = result + str(nELE)
							result = result + ')'
							result = result + eleHeavy
							result = result + '('
							result = result + str(nELE)
							result = result + ')'

			return result

	def getDictComposition(self, inputComposition):

		outputDictComposition = {}  # 这是一个以名字为key，以个数为number的dict

		LastIndexLeftBracket = 0
		LastIndexRightBracket = -1  # 第一开始，+1后从0开始

		lengthInputString = len(inputComposition)
		for i in range(lengthInputString):

			if inputComposition[i] == '(':
				eleName = inputComposition[LastIndexRightBracket + 1:i]
				LastIndexLeftBracket = i

			if inputComposition[i] == ')':

				eleAtomNum = int(inputComposition[LastIndexLeftBracket + 1:i])

				if outputDictComposition.__contains__(eleName):

					outputDictComposition[eleName] = outputDictComposition[eleName] + eleAtomNum  # eleAtomNum可能是负数

				else:

					outputDictComposition[eleName] = eleAtomNum

				LastIndexRightBracket = i

		return outputDictComposition

	def getStrComposition(self, inputSeq, inputMod, inputGLC, inputLIK, inputINI):

		result = 'H(2)O(1)'  # 写死的

		for aa in inputSeq:

			if inputINI.DICT1_AA_COM.__contains__(aa):

				result = result + inputINI.DICT1_AA_COM[aa]

			else:

				continue

		# modification
		nMOD = toolCountCharInString(inputMod, ';')  # 标准形式为6,Carbamidomethyl[C];14,Carbamidomethyl[C];
		try:
			for iMOD in range(nMOD):
				nameMOD = toolGetWord(toolGetWord(inputMod, iMOD, ';'), 1, ',')
				result = result + inputINI.DICT2_MOD_COM[nameMOD]
		except:
			logToUser(INFO_TO_USER_FunctionComposition[0] + inputMod)
			logGetError(INFO_TO_USER_FunctionComposition[1])

		# glyco
		nGLC = toolCountCharInString(inputGLC, ';')
		try:
			for iGLC in range(nGLC):
				nameGLC = toolGetWord(toolGetWord(inputGLC, iGLC, ';'), 1, ',')
				if inputINI.DICT3_GLYCO_COM.__contains__(nameGLC):
					result = result + inputINI.DICT3_GLYCO_COM[nameGLC]
				else:
					logGetWarning(inputGLC)
		except:
			logToUser(INFO_TO_USER_FunctionComposition[0] + inputGLC)
			logGetError(INFO_TO_USER_FunctionComposition[2])

		# linker
		nLINK = toolCountCharInString(inputLIK, ';')  # 标准形式为同修饰
		try:
			for iLINK in range(nLINK):
				nameLINK = toolGetWord(toolGetWord(inputLIK, iLINK, ';'), 1, ',')
				result = result + inputINI.DICT4_LINKER_COM[nameLINK]
			if inputSeq.find('-') != -1:
				result += 'H(2)O(1)'
		except:
			logToUser(INFO_TO_USER_FunctionComposition[0] + inputLIK)
			logGetError(INFO_TO_USER_FunctionComposition[3])

		return result

