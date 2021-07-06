# encoder stolen from phil

# import stuff

def hotEncodingAAString(myString):
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	char_to_int = dict((c, i) for i, c in enumerate(alphabet))
	int_to_char = dict((i, c) for i, c in enumerate(alphabet))
	onehot_encoded = list()
	integer_encoded = [char_to_int[char] for char in myString]
	for value in integer_encoded:
		letter = [0 for _ in range(len(alphabet))]
		letter[value] = 1
		onehot_encoded.append(letter)
	#print(onehot_encoded)
	return[onehot_encoded]


def hotEncodingAAStringflat(myString):
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	char_to_int = dict((c, i) for i, c in enumerate(alphabet))
	int_to_char = dict((i, c) for i, c in enumerate(alphabet))
	onehot_encoded = list()
	integer_encoded = [char_to_int[char] for char in myString]
	for value in integer_encoded:
		letter = [0 for _ in range(len(alphabet))]
		letter[value] = 1
		onehot_encoded.append(letter)
	#print(onehot_encoded)
	onehot_encoded = sum(onehot_encoded, [])
	return[onehot_encoded]



def batchhotEncodingAAStringflat(myStringList):
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	char_to_int = dict((c, i) for i, c in enumerate(alphabet))
	int_to_char = dict((i, c) for i, c in enumerate(alphabet))
	onehot_encodeds = []
	for myString in myStringList:
		onehot_encoded = list()
		integer_encoded = [char_to_int[char] for char in myString]
		for value in integer_encoded:
			letter = [0 for _ in range(len(alphabet))]
			letter[value] = 1
			onehot_encoded.append(letter)
		#print(onehot_encoded)
		onehot_encoded = sum(onehot_encoded, [])
		onehot_encodeds.append(onehot_encoded)
	return onehot_encodeds


def label_binarizer(listoflabels):
	'''
	transform labels (Binder, NonBinder) to binary labels (1,1)
	:param listoflabels:
	:return:
	'''
	binaries = []
	for label in listoflabels:
		if label == 'Binder':
			binary = 1
		else:
			binary = 0
		binaries.append(binary)
	return binaries

# run stuff
a  = hotEncodingAAString('VICTR')
b  = batchhotEncodingAAStringflat(['VICTR', 'DLA'])

print(a, len(a))
print(b, len(b))
bins = label_binarizer(['Binder', 'NonBinder'])
print(bin)